#!/bin/bash

# ========================================
# Pipeline automatisé pour l'appel de variants :)
# Auteur: [Votre nom]
# Date: [Date actuelle]
# ========================================

#Chargement des modules permettant d'exécuter les commandes du script
module load sabre fastqc cutadapt parallel bwa samtools bcftools snpEff

# ==== Vérifier les arguments de ligne de commande ====
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fastq_file> <barcode_file>  <genome_ref_file>"
    exit 1
fi

# ==== Définir les variables ====
DATA=$1
BARCODE=$2
TOOL="scripts/sabre.sh"
ADAP="AGATCGGAA"
REF=$3
CPU=4
OUT_FASTQC_PRE="results/fastqc_pre_trim"
OUT_FASTQC_POST="results/fastqc_post_trim"
OUT_ASSEMBLY="results/assembly"
OUT_RESULTS="results/results1"

# ==== Vérifications initiales ====
if [ ! -f "$DATA" ]; then
    echo "ERREUR : Fichier FASTQ $DATA introuvable." >&2
    exit 1
fi

if [ ! -f "$BARCODE" ]; then
    echo "ERREUR : Fichier de barcodes $BARCODE introuvable." >&2
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "ERREUR : Génome de référence $REF introuvable." >&2
    exit 1
fi

# ==== Étape 1 : Démultiplexage ====
exec &> logs/sabre.log
mkdir -p results
echo "Démarrage du démultiplexage avec sabre..."
$TOOL se -f $DATA -b $BARCODE -u results/unk.fastq
if [ $? -ne 0 ]; then
    echo "Échec de l'étape de démultiplexage avec sabre." >&2
    exit 1
fi
echo "Démultiplexage terminé avec succès."

# ==== Étape 2 : Analyse de qualité pré-trim ====
mkdir -p $OUT_FASTQC_PRE
echo "Analyse de qualité (pré-trim) avec FastQC..."
parallel -j $CPU fastqc -o $OUT_FASTQC_PRE {} ::: $DATA
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'analyse de qualité pré-trim." >&2
    exit 1
fi
echo "Analyse de qualité pré-trim terminée."

# ==== Étape 3 : Trimming ====
exec &> logs/cutadapt.log
echo "Démarrage du trimming avec Cutadapt..."
mkdir -p trimmed
parallel -j $CPU cutadapt -a $ADAP -m 50 -o trimmed/{}.fastq {} ::: $(ls -1 data/*.fq.gz)
if [ $? -ne 0 ]; then
    echo "Erreur lors du trimming." >&2
    exit 1
fi
echo "Trimming terminé."

# ==== Étape 4 : Analyse de qualité post-trim ====
mkdir -p $OUT_FASTQC_POST
echo "Analyse de qualité (post-trim) avec FastQC..."
parallel -j $CPU fastqc -o $OUT_FASTQC_POST {} ::: trimmed/*.fastq
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'analyse de qualité post-trim." >&2
    exit 1
fi
echo "Analyse de qualité post-trim terminée."

# ==== Étape 5 : Alignement ====
exec &> logs/bwa.log
mkdir -p aligned
echo "Alignement avec BWA..."
parallel -j $CPU bwa mem -t $CPU $REF {} ">" aligned/{}.sam ::: trimmed/*.fastq
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'alignement." >&2
    exit 1
fi
echo "Alignement terminé."

# ==== Étape 6 : Conversion SAM -> BAM ====
exec &> logs/samtools.log
echo "Conversion des fichiers SAM en BAM..."
mkdir -p bam_files
parallel -j $CPU samtools view -bS aligned/{}.sam ">" bam_files/{}.bam ::: $(ls -1 aligned/*.sam | sed 's/\.sam//')
parallel -j $CPU samtools sort bam_files/{}.bam -o bam_files/{}.sorted.bam ::: $(ls -1 bam_files/*.bam | sed 's/\.bam//')
parallel -j $CPU samtools index bam_files/{}.sorted.bam ::: $(ls -1 bam_files/*.sorted.bam)
echo "Conversion SAM -> BAM terminée."

# ==== Étape 7 : Appel de variants ====
exec &> logs/variant_calling.log
mkdir -p $OUT_RESULTS
echo "Appel de variants avec SAMtools et BCFtools..."
samtools mpileup -g -f $REF -b bam_files/*.sorted.bam > $OUT_RESULTS/variants.bcf
bcftools call -mv -O v -o $OUT_RESULTS/variants.vcf $OUT_RESULTS/variants.bcf
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'appel de variants." >&2
    exit 1
fi
echo "Appel de variants terminé."

# ==== Étape 8 (Bonus) : Annotation ====
exec &> logs/annotation.log
echo "Annotation des variants (optionnelle)..."
if command -v snpEff &> /dev/null; then
    snpEff annotate -v $OUT_RESULTS/variants.vcf > $OUT_RESULTS/variants_annotated.vcf
    echo "Annotation terminée avec succès."
else
    echo "SnpEff non installé, annotation non effectuée."
fi

echo "Pipeline terminé avec succès."