#!/bin/bash
# Verbose Genomics Pipeline with Detailed Logging
# Exit immediately on any error
set -eE

# Trap any errors and provide detailed error information
trap 'handle_error "$BASH_COMMAND" $?' ERR

# Comprehensive error handling function
handle_error() {
    local command="$1"
    local exit_code="$2"
    
    echo "❌ ERROR OCCURRED ❌"
    echo "Command: $command"
    echo "Exit Code: $exit_code"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    
    # Log detailed error information
    {
        echo "------- ERROR DETAILS -------"
        echo "Command: $command"
        echo "Exit Code: $exit_code"
        echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Current Directory: $(pwd)"
        echo "Script Arguments: $@"
    } >> "$OUTPUT_DIR/pipeline_error.log"
    
    # Optional: Send an email or notification about the error
    # Uncomment and configure if needed
    # send_error_notification "$command" "$exit_code"
    
    exit 1
}

# Determine script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configurable paths with environment variable support
DATA_DIR="${DATA_DIR:-$SCRIPT_DIR/data}"
REFERENCE_DIR="${REFERENCE_DIR:-$SCRIPT_DIR/reference}"
OUTPUT_DIR="${OUTPUT_DIR:-$SCRIPT_DIR/results}"
TOOLS_DIR="${TOOLS_DIR:-$SCRIPT_DIR/tools}"

# Verbose logging function
verbose_log() {
    local message="$1"
    local log_file="${OUTPUT_DIR}/pipeline_verbose.log"
    
    # Create log directory if it doesn't exist
    mkdir -p "$(dirname "$log_file")"
    
    # Print to console with color
    echo -e "\e[34m[PIPELINE]\e[0m \e[32m$message\e[0m"
    
    # Log to file
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$log_file"
}

# Progress tracking function
track_progress() {
    local stage="$1"
    verbose_log "🔹 STAGE: $stage"
}

# Tool path detection with explicit error handling
detect_tool() {
    local tool_name="$1"
    local tool_path
    
    tool_path=$(which "$tool_name" 2>/dev/null)
    
    if [ -z "$tool_path" ]; then
        echo "❌ ERROR: $tool_name not found in system PATH"
        exit 1
    fi
    
    echo "$tool_path"
}

# Detect and set tool paths
FASTQC_PATH=$(detect_tool fastqc)
CUTADAPT_PATH=$(detect_tool cutadapt)
BWA_PATH=$(detect_tool bwa)
SAMTOOLS_PATH=$(detect_tool samtools)
BCFTOOLS_PATH=$(detect_tool bcftools)
JAVA_PATH=$(detect_tool java)

# Genome Configuration
GENOME_NAME="Gmax_275"
REFERENCE_GENOME="${REFERENCE_DIR}/${GENOME_NAME}/reference_genome.fasta"
SNPEFF_PATH="${TOOLS_DIR}/snpEff"
SNPEFF_DB="${SNPEFF_PATH}/data/Glycine_max"

# Input Validation Function
validate_inputs() {
    track_progress "Input Validation"
    
    # Check each input file
    local input_fastq="$1"
    local barcode_file="$2"
    local reference_genome="$3"
    
    for file in "$input_fastq" "$barcode_file" "$reference_genome"; do
        if [ ! -f "$file" ]; then
            echo "❌ MISSING FILE: $file"
            exit 1
        fi
    done
    
    verbose_log "✅ All input files validated successfully"
}

# Function to prepare BWA index
prepare_bwa_index() {
    local reference_genome="$1"
    local base_genome="${reference_genome%.*}"
    
    # Check if index files exist
    if [ ! -f "${reference_genome}.bwt" ] && [ ! -f "${base_genome}.bwt" ]; then
        track_progress "Generating BWA Index Files"
        verbose_log "🔧 Creating BWA index for reference genome"
        
        # Try to create index
        if ! bwa index "$reference_genome"; then
            echo "❌ ERROR: Failed to create BWA index for $reference_genome"
            exit 1
        fi
        
        verbose_log "✅ BWA index created successfully"
    else
        verbose_log "✅ BWA index files already exist"
    fi
}

# Main Pipeline Function
run_genomics_pipeline() {
    local input_fastq="$1"
    local barcode_file="$2"
    local reference_genome="$3"
    
    # Validate inputs first
    validate_inputs "$input_fastq" "$barcode_file" "$reference_genome"
    
    # Create output directories
    track_progress "Creating Output Directories"
    mkdir -p "$OUTPUT_DIR"/{fastqc_pre_trim,fastqc_post_trim,trimmed_reads,alignment,variants,snpeff}
    
    # Pre-Trimming Quality Control
    track_progress "Pre-Trimming Quality Control"
    "$FASTQC_PATH" -t 4 -o "$OUTPUT_DIR/fastqc_pre_trim" "$input_fastq"
    verbose_log "✅ Pre-trimming FastQC complete"
    
    # Adapter Trimming
    track_progress "Adapter Trimming"
    "$CUTADAPT_PATH" -a AGATCGGAA -m 50 -o "$OUTPUT_DIR/trimmed_reads/trimmed.fastq" "$input_fastq"
    verbose_log "✅ Reads trimmed successfully"
    
    # Post-Trimming Quality Control
    track_progress "Post-Trimming Quality Control"
    "$FASTQC_PATH" -t 4 -o "$OUTPUT_DIR/fastqc_post_trim" "$OUTPUT_DIR/trimmed_reads/trimmed.fastq"
    verbose_log "✅ Post-trimming FastQC complete"
    prepare_bwa_index "$reference_genome"    

    # Alignment
    track_progress "Genome Alignment"
    "$BWA_PATH" mem -t 4 "$reference_genome" "$OUTPUT_DIR/trimmed_reads/trimmed.fastq" > "$OUTPUT_DIR/alignment/aligned.sam"
    verbose_log "✅ Alignment completed"
    
    # SAM to BAM Conversion
    track_progress "SAM to BAM Conversion"
    "$SAMTOOLS_PATH" view -b "$OUTPUT_DIR/alignment/aligned.sam" > "$OUTPUT_DIR/alignment/aligned.bam"
    "$SAMTOOLS_PATH" sort "$OUTPUT_DIR/alignment/aligned.bam" -o "$OUTPUT_DIR/alignment/sorted.bam"
    "$SAMTOOLS_PATH" index "$OUTPUT_DIR/alignment/sorted.bam"
    verbose_log "✅ BAM files processed successfully"
    
    # Variant Calling
    track_progress "Variant Calling"
    "$BCFTOOLS_PATH" mpileup -f "$reference_genome" "$OUTPUT_DIR/alignment/sorted.bam" > "$OUTPUT_DIR/variants/variants.bcf"
    "$BCFTOOLS_PATH" call -mv "$OUTPUT_DIR/variants/variants.bcf" > "$OUTPUT_DIR/variants/variants.vcf"
    verbose_log "✅ Variant calling completed"
    
    # Variant Annotation
    track_progress "Variant Annotation"
    "$JAVA_PATH" -jar "$SNPEFF_PATH/snpEff.jar" -v -o vcf "$GENOME_NAME" "$OUTPUT_DIR/variants/variants.vcf" > "$OUTPUT_DIR/snpeff/annotated_variants.vcf"
    verbose_log "✅ Variant annotation finished"
    
    echo "🎉 Pipeline completed successfully 🎉"
}

# Script Execution
main() {
    # Check for correct number of arguments
    if [ $# -ne 3 ]; then
        echo "Usage: $0 <input_fastq> <barcode_file> <reference_genome>"
        exit 1
    fi
    
    # Start the pipeline
    run_genomics_pipeline "$@"
}

# Run main with all script arguments
main "$@"
