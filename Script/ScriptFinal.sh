#!/bin/bash
# ========================================
# Universal Pipeline for Variant Calling
# ========================================

# Load required modules
module load sabre
module load fastqc
module load cutadapt
module load bwa
module load samtools
module load bcftools

# Explicitly set paths for directories and files
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
DATA_DIR="$SCRIPT_DIR/../Data"
LOG_DIR="$SCRIPT_DIR/../Logs"
echo "DEBUG: SCRIPT_DIR is $SCRIPT_DIR"
echo "DEBUG: DATA_DIR is $DATA_DIR"
echo "DEBUG: LOG_DIR is $LOG_DIR"

# Create a Logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Check if Data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "ERREUR: Répertoire de données introuvable: $DATA_DIR" | tee -a "$LOG_DIR/errors.log"
    exit 1
fi

# Locate required files
INPUT_FASTQ=$(find "$DATA_DIR" -type f -name "*.fq.gz" | head -n 1)
BARCODE_FILE=$(find "$DATA_DIR" -type f -name "*.txt" | head -n 1)
REF_GENOME=$(find "$DATA_DIR" -type f -name "*.fa.gz" | head -n 1)

# Ensure required files are found
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

# Debugging information
echo "DEBUG: INPUT_FASTQ = $INPUT_FASTQ" | tee -a "$LOG_DIR/debug.log"
echo "DEBUG: BARCODE_FILE = $BARCODE_FILE" | tee -a "$LOG_DIR/debug.log"
echo "DEBUG: REF_GENOME = $REF_GENOME" | tee -a "$LOG_DIR/debug.log"

# Variables
CPU=4
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
BASE_DIR="$SCRIPT_DIR/../Résultats/GBS_data_${TIMESTAMP}"

# Create output directories
mkdir -p "$BASE_DIR/fastqc_pre_trim"
mkdir -p "$BASE_DIR/fastqc_post_trim"
mkdir -p "$BASE_DIR/trimmed_reads"
mkdir -p "$BASE_DIR/alignment"
mkdir -p "$BASE_DIR/results"
mkdir -p "$BASE_DIR/logs"

# Log and error-handling function
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

# Logs
exec > >(tee -a "$LOG_DIR/pipeline.log") 2>&1

# Decompress the reference genome if necessary
if [[ "$REF_GENOME" == *.gz ]]; then
    log_and_check "gunzip -c $REF_GENOME | bgzip > ${REF_GENOME%.gz}.bgz" \
        "Decompression of reference genome failed"
    REF_GENOME="${REF_GENOME%.gz}.bgz"
fi

# Step 1: Index Reference Genome (if needed)
if [ ! -f "$REF_GENOME.bwt" ]; then
    echo "Indexing reference genome with BWA..."
    log_and_check "bwa index $REF_GENOME" "Échec de l'indexation du génome de référence."
    echo "Indexing completed successfully."
fi

# Step 2: Demultiplexing with Sabre
log_and_check "sabre se -f $INPUT_FASTQ -b $BARCODE_FILE -u $BASE_DIR/unknown.fastq" \
    "Échec du démultiplexage avec sabre"

# Step 3: FastQC Pre-trim
log_and_check "fastqc -o $BASE_DIR/fastqc_pre_trim $INPUT_FASTQ" \
    "Échec de l'analyse FastQC pré-trim"

# Step 4: Trimming with Cutadapt
TRIMMED_FASTQ="$BASE_DIR/trimmed_reads/$(basename "$INPUT_FASTQ" .fq.gz)_trimmed.fq.gz"
log_and_check "cutadapt -a AGATCGGAA -m 50 -o $TRIMMED_FASTQ $INPUT_FASTQ" \
    "Échec du trimming avec Cutadapt"

# Step 5: FastQC Post-trim
log_and_check "fastqc -o $BASE_DIR/fastqc_post_trim $TRIMMED_FASTQ" \
    "Échec de l'analyse FastQC post-trim"

# Step 6: Alignment with BWA
SAM_FILE="$BASE_DIR/alignment/aligned.sam"
log_and_check "bwa mem -t $CPU $REF_GENOME $TRIMMED_FASTQ > $SAM_FILE" \
    "Échec de l'alignement BWA"

# Step 7: SAM to BAM Conversion with Samtools
BAM_FILE="$BASE_DIR/alignment/aligned.bam"
SORTED_BAM="$BASE_DIR/alignment/sorted.bam"
log_and_check "samtools view -b $SAM_FILE > $BAM_FILE" \
    "Échec de la conversion SAM à BAM"

log_and_check "samtools sort $BAM_FILE -o $SORTED_BAM" \
    "Échec du tri BAM"

log_and_check "samtools index $SORTED_BAM" \
    "Échec de l'indexation BAM"

# Step 8: Variant Calling with BCFtools  
VARIANTS_BCF="$BASE_DIR/results/variants.bcf"  
VARIANTS_VCF="$BASE_DIR/results/variants.vcf"  

# Use bcftools mpileup instead of samtools mpileup  
log_and_check "bcftools mpileup -Ou -f $REF_GENOME $SORTED_BAM | bcftools call -mv -Ov -o $VARIANTS_VCF" \
    "Échec du mpileup et de l'appel de variants"

echo "Pipeline terminé avec succès. Résultats dans $BASE_DIR" | tee -a "$LOG_DIR/pipeline.log"
