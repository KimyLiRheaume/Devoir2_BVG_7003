#!/bin/bash
#!/bin/bash

# Define relative paths
DATA_DIR="./Data"
RESULTS_DIR="./Résultats"
SCRIPT_DIR="./Script"
LOGS_DIR="./Logs"
README_FILE="./README.md"

# Example operations

# List files in Data directory
echo "Listing files in Data directory:"
ls $DATA_DIR

# List files in Results directory
echo "Listing files in Resultats directory:"
ls $RESULTS_DIR

# Run a script from the Script directory
echo "Running a script from the Script directory:"
bash $SCRIPT_DIR/Script_v4.sh

# Create a log file in Logs directory
echo "Creating a log file in Logs directory:"
echo "Log entry" > $LOGS_DIR/your_log_file.log

# Display the content of README.md
#echo "Displaying the content of README.md:"
#cat $README_FILE

# Example of processing a file from Data directory and saving the result in Results directory
echo "Processing a file from Data directory and saving the result in Results directory:"
cp $DATA_DIR/your_file.fastq $RESULTS_DIR/processed_file.fastq

echo "Script execution completed."


# Check if required arguments are provided
if [ $# -ne 3 ]; then
    echo "Utilisation: $0 <fichier_fastq> <fichier_barcode> <reference_genome>"
    exit 1
fi

# Arguments
INPUT_FASTQ="$1"
BARCODE_FILE="$2"
REF_GENOME="$3"

# Définir les variables
CPU=4
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
BASE_DIR="GBS_data_${TIMESTAMP}"
SABRE_PATH="../sabre/sabre" #Possiblement a changer pour le path sur votre ordi

# Créer des répertoires de sortie
mkdir -p "$BASE_DIR"
mkdir -p "$BASE_DIR/fastqc_pre_trim"
mkdir -p "$BASE_DIR/fastqc_post_trim"
mkdir -p "$BASE_DIR/trimmed_reads"
mkdir -p "$BASE_DIR/alignment"
mkdir -p "$BASE_DIR/results"

# Fonction de log et de gestion des erreurs
log_and_check() {
    local command="$1"
    local error_message="$2"
    
    echo "Exécution: $command"
    eval "$command"
    
    if [ $? -ne 0 ]; then
        echo "ERREUR: $error_message"
        exit 1
    fi
}

# Logs
exec > >(tee "$BASE_DIR/pipeline.log") 2>&1

# Vérifier l'existence des fichiers d'entrée
for file in "$INPUT_FASTQ" "$BARCODE_FILE" "$REF_GENOME"; do
    if [ ! -f "$file" ]; then
        echo "ERREUR: Fichier manquant - $file"
        exit 1
    fi
done

# Vérifier l'existence de sabre
if [ ! -x "$SABRE_PATH" ]; then
    echo "ERREUR: Executable sabre non trouvé à $SABRE_PATH"
    exit 1
fi

# Démultiplexage 
log_and_check "$SABRE_PATH se -f $INPUT_FASTQ -b $BARCODE_FILE -u $BASE_DIR/unknown.fastq" \
    "Échec du démultiplexage"

# Décompresser le génome de référence si nécessaire
if [[ "$REF_GENOME" == *.gz ]]; then
    log_and_check "gunzip -k $REF_GENOME" \
        "Échec de la décompression du génome de référence"
    REF_GENOME="${REF_GENOME%.gz}"
fi

# Indexer le génome de référence si l'index n'existe pas
if [ ! -f "${REF_GENOME}.bwt" ]; then
    log_and_check "bwa index $REF_GENOME" \
        "Échec de l'indexation du génome de référence"
fi

# FastQC pré-trimming
log_and_check "fastqc -o $BASE_DIR/fastqc_pre_trim $INPUT_FASTQ" \
    "Échec de FastQC pré-trimming"

# Trimming avec cutadapt
log_and_check "cutadapt -a AGATCGGAA -m 50 -o $BASE_DIR/trimmed_reads/trimmed.fastq $INPUT_FASTQ" \
    "Échec du trimming avec cutadapt"

# FastQC post-trimming
log_and_check "fastqc -o $BASE_DIR/fastqc_post_trim $BASE_DIR/trimmed_reads/trimmed.fastq" \
    "Échec de FastQC post-trimming"

# Alignement BWA
log_and_check "bwa mem -t $CPU $REF_GENOME $BASE_DIR/trimmed_reads/trimmed.fastq > $BASE_DIR/alignment/aligned.sam" \
    "Échec de l'alignement BWA"

# Conversion SAM à BAM
log_and_check "samtools view -b $BASE_DIR/alignment/aligned.sam > $BASE_DIR/alignment/aligned.bam" \
    "Échec de la conversion SAM à BAM"

log_and_check "samtools sort $BASE_DIR/alignment/aligned.bam -o $BASE_DIR/alignment/sorted.bam" \
    "Échec du tri BAM"

log_and_check "samtools index $BASE_DIR/alignment/sorted.bam" \
    "Échec de l'indexation BAM"

# Indexer le génome de référence pour samtools
if [ ! -f "$REF_GENOME.fai" ]; then
    log_and_check "samtools faidx $REF_GENOME" \
        "Échec de l'indexation du génome pour samtools"
fi

# Appel de variants
log_and_check "samtools mpileup -g -f $REF_GENOME $BASE_DIR/alignment/sorted.bam > $BASE_DIR/results/variants.bcf" \
    "Échec de l'étape mpileup"

log_and_check "bcftools call -mv $BASE_DIR/results/variants.bcf > $BASE_DIR/results/variants.vcf" \
    "Échec de l'appel de variants"

echo "Processus terminé avec succès dans $BASE_DIR"
