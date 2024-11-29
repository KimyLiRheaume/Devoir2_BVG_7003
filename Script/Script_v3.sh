#!/bin/bash

# ========================================
# Pipeline automatisé pour l'appel de variants
# Auteur: [Votre nom]
# Date: [Date actuelle]
# ========================================

# Pipeline automatisé pour l'appel de variants
# Auteur: [Votre nom]
# Date: [Date actuelle]
# ========================================

# Directory for SnpEff
mkdir -p ~/Devoir2_BVG_7003/Script/snpEff && cd ~/Devoir2_BVG_7003/Scripts/snpEff

# Download the latest version of SnpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Extract the downloaded zip file
unzip snpEff_latest_core.zip

# Add SnpEff to PATH
echo 'export PATH=$HOME/Devoir2_BVG_7003/Script/snpEff:$PATH' >> ~/.bashrc
source ~/.bashrc

# Test the installation
java -jar ~/Devoir2_BVG_7003/Script/snpEff/snpEff.jar -h

# (Optional) Download a genome database
java -jar ~/Devoir2_BVG_7003/Script/snpEff/snpEff.jar download GRCh38.99

#Directory for Parallel
mkdir -p ~/Devoir2_BVG_7003/Script/paralell && cd ~/Devoir2_BVG_7003/Scripts/paralell

# Download the latest version of Parallel
(wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
pi.dk

# Add paralell to PATH
echo 'export PATH=$HOME/Devoir2_BVG_7003/Script/parallel:$PATH' >> ~/.bashrc
source ~/.bashrc

#module load python/3.7 sabre fastqc cutadapt parallel bwa samtools bcftools snpEff

# ==== Vérifier les arguments de ligne de commande ====
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fastq_file> <barcode_file>"
    exit 1
fi

# ==== Définir les variables ====
DATA=$1
BARCODE=$2
TOOL="sabre"
ADAP="AGATCGGAA"
REF="Data/Gmax_275_v2.0.fa.gz"
CPU=4
OUT_FASTQC_PRE="Résultats/fastqc_pre_trim"
OUT_FASTQC_POST="Résultats/fastqc_post_trim"
OUT_ASSEMBLY="Résultats/assembly"
OUT_RESULTS="Résultats/results1"

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
mkdir -p Résultats
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
