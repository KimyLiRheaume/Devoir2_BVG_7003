#!/bin/bash
# ========================================
# Pipeline Universel pour l'Appel de Variants
# ========================================

# Charger les modules nécessaires
module load sabre
module load fastqc
module load cutadapt
module load bwa
module load samtools
module load bcftools

# Définir explicitement les chemins des répertoires et fichiers
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
DATA_DIR="$SCRIPT_DIR/../Data"
LOG_DIR="$SCRIPT_DIR/../Logs"
echo "DEBUG: SCRIPT_DIR is $SCRIPT_DIR"
echo "DEBUG: DATA_DIR is $DATA_DIR"
echo "DEBUG: LOG_DIR is $LOG_DIR"

# Créer un répertoire Logs s'il n'existe pas
mkdir -p "$LOG_DIR"

# Vérifier si le répertoire Data existe
if [ ! -d "$DATA_DIR" ]; then
    echo "ERREUR: Répertoire de données introuvable: $DATA_DIR" | tee -a "$LOG_DIR/errors.log"
    exit 1
fi

# Localiser les fichiers nécessaires
INPUT_FASTQ=$(find "$DATA_DIR" -type f -name "*.fq.gz" | head -n 1)
BARCODE_FILE=$(find "$DATA_DIR" -type f -name "*.txt" | head -n 1)
REF_GENOME=$(find "$DATA_DIR" -type f -name "*.fa.gz" | head -n 1)

# Vérifier que les fichiers requis sont trouvés
if [ -z "$INPUT_FASTQ" ]; then
    echo "ERREUR: Aucun fichier FASTQ (*.fq.gz) trouvé dans $DATA_DIR." | tee -a "$LOG_DIR/errors.log"
    exit 1
fi
if [ -z "$BARCODE_FILE" ]; then
    echo "ERREUR: Aucun fichier de barcodes (*.txt) trouvé dans $DATA_DIR." | tee -a "$LOG_DIR/errors.log"
    exit 1
fi
if [ -z "$REF_GENOME" ]; then
    echo "ERREUR: Aucun fichier de génome de référence (*.fa.gz) trouvé dans $DATA_DIR." | tee -a "$LOG_DIR/errors.log"
    exit 1
fi

# Informations de débogage
echo "DEBUG: INPUT_FASTQ = $INPUT_FASTQ" | tee -a "$LOG_DIR/debug.log"
echo "DEBUG: BARCODE_FILE = $BARCODE_FILE" | tee -a "$LOG_DIR/debug.log"
echo "DEBUG: REF_GENOME = $REF_GENOME" | tee -a "$LOG_DIR/debug.log"

# Variables
CPU=4
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
BASE_DIR="$SCRIPT_DIR/../Résultats/GBS_data_${TIMESTAMP}"

# Créer les répertoires de sortie
mkdir -p "$BASE_DIR/fastqc_pre_trim"
mkdir -p "$BASE_DIR/fastqc_post_trim"
mkdir -p "$BASE_DIR/trimmed_reads"
mkdir -p "$BASE_DIR/alignment"
mkdir -p "$BASE_DIR/results"
mkdir -p "$BASE_DIR/logs"

# Fonction de journalisation et de gestion des erreurs
log_and_check() {
    local command="$1"
    local error_message="$2"
    local log_file="$LOG_DIR/pipeline.log"

    echo "Exécution: $command" | tee -a "$log_file"
    eval "$command" &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "ERREUR: $error_message" | tee -a "$LOG_DIR/errors.log"
        exit 1
    fi
}

# Journaux
exec > >(tee -a "$LOG_DIR/pipeline.log") 2>&1

# Décompresser le génome de référence si nécessaire
if [[ "$REF_GENOME" == *.gz ]]; then
    log_and_check "gunzip -c $REF_GENOME | bgzip > ${REF_GENOME%.gz}.bgz" \
        "Échec de la décompression du génome de référence"
    REF_GENOME="${REF_GENOME%.gz}.bgz"
fi

# Étape 1 : Indexation du génome de référence
if [ ! -f "$REF_GENOME.bwt" ]; then
    echo "Indexation du génome de référence avec BWA..."
    log_and_check "bwa index $REF_GENOME" "Échec de l'indexation du génome de référence."
    echo "Indexation terminée avec succès."
fi

# Étape 2 : Démultiplexage avec Sabre
log_and_check "sabre se -f $INPUT_FASTQ -b $BARCODE_FILE -u $BASE_DIR/unknown.fastq" \
    "Échec du démultiplexage avec sabre"

# Étape 3 : Analyse FastQC avant trimming
log_and_check "fastqc -o $BASE_DIR/fastqc_pre_trim $INPUT_FASTQ" \
    "Échec de l'analyse FastQC avant trimming"

# Étape 4 : Trimming avec Cutadapt
TRIMMED_FASTQ="$BASE_DIR/trimmed_reads/$(basename "$INPUT_FASTQ" .fq.gz)_trimmed.fq.gz"
log_and_check "cutadapt -a AGATCGGAA -m 50 -o $TRIMMED_FASTQ $INPUT_FASTQ" \
    "Échec du trimming avec Cutadapt"

# Étape 5 : Analyse FastQC après trimming
log_and_check "fastqc -o $BASE_DIR/fastqc_post_trim $TRIMMED_FASTQ" \
    "Échec de l'analyse FastQC après trimming"

# Étape 6 : Alignement avec BWA
SAM_FILE="$BASE_DIR/alignment/aligned.sam"
log_and_check "bwa mem -t $CPU $REF_GENOME $TRIMMED_FASTQ > $SAM_FILE" \
    "Échec de l'alignement avec BWA"

# Étape 7 : Conversion SAM en BAM avec Samtools
BAM_FILE="$BASE_DIR/alignment/aligned.bam"
SORTED_BAM="$BASE_DIR/alignment/sorted.bam"
log_and_check "samtools view -b $SAM_FILE > $BAM_FILE" \
    "Échec de la conversion SAM en BAM"

log_and_check "samtools sort $BAM_FILE -o $SORTED_BAM" \
    "Échec du tri BAM"

log_and_check "samtools index $SORTED_BAM" \
    "Échec de l'indexation BAM"

# Étape 8 : Appel de variants avec BCFtools
VARIANTS_BCF="$BASE_DIR/results/variants.bcf"
VARIANTS_VCF="$BASE_DIR/results/variants.vcf"

log_and_check "bcftools mpileup -Ou -f $REF_GENOME $SORTED_BAM | bcftools call -mv -Ov -o $VARIANTS_VCF" \
    "Échec du mpileup et de l'appel de variants"

echo "Pipeline terminé avec succès. Résultats dans $BASE_DIR" | tee -a "$LOG_DIR/pipeline.log"
