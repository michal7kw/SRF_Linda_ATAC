#!/bin/bash
# =============================================================================
# SCRIPT: 01_extract_barcodes.sh
# PURPOSE: Extract cell barcodes from R2 reads for scATAC-seq analysis
# 
# FUNCTIONALITY:
# - Extracts rightmost 16bp from 24bp R2 reads (positions 9-24)
# - Handles 10x Genomics barcode structure: [8bp spacer][16bp barcode]
# - Processes multiple samples using SLURM array jobs
# - Validates extracted barcode length and count
# 
# INPUT: Raw FASTQ files (R2 reads containing barcodes)
# OUTPUT: Processed FASTQ files with extracted 16bp barcodes
# 
# EDGE CASES:
# - Skips extraction if output file already exists
# - Validates barcode length to ensure correct extraction
# - Handles gzipped input/output files
# =============================================================================

#SBATCH --job-name=extract_barcodes
#SBATCH --output=logs/01_extract_barcodes_%a.out
#SBATCH --error=logs/01_extract_barcodes_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq

# =============================================================================
# ENVIRONMENT SETUP
# PURPOSE: Initialize conda environment and set strict error handling
# 
# OPERATIONS:
# - Activates conda environment with required bioinformatics tools
# - Enables strict bash error handling (exit on error, undefined vars, pipe failures)
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# CONFIGURATION AND SAMPLE SELECTION
# PURPOSE: Define sample names and directory paths for processing
# 
# OPERATIONS:
# - Maps SLURM array task ID to specific sample name
# - Sets input directory containing raw FASTQ files
# - Sets output directory for processed barcode files
# 
# INPUT EXPECTATIONS:
# - SLURM_ARRAY_TASK_ID: 0 or 1 (corresponding to control/mutant samples)
# - Raw FASTQ files in DATA_DIR with naming pattern: {SAMPLE}_R2_001.fastq.gz
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/ATAC_data/DATA/nestin_processed"

echo "========================================="
echo "Step 1: Extracting barcodes for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY CREATION
# PURPOSE: Ensure output directory exists for processed files
# 
# OPERATIONS:
# - Creates directory structure if it doesn't exist
# - Uses -p flag to create parent directories and avoid errors if exists
# =============================================================================
mkdir -p "$PROCESSED_DATA_DIR"

# =============================================================================
# BARCODE EXTRACTION SETUP
# PURPOSE: Define output file and explain 10x Genomics barcode structure
# 
# ALGORITHM:
# - 10x Genomics R2 structure: [8bp spacer][16bp cell barcode] = 24bp total
# - Extract rightmost 16bp (positions 9-24) containing actual cell barcode
# - Spacer region (positions 1-8) is discarded as it's not part of barcode
# 
# OUTPUT: FASTQ file with extracted 16bp barcodes maintaining read structure
# =============================================================================
EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

echo "DEBUG: R2 structure according to facility:"
echo "  [8bp SPACER][16bp 10x barcode] = 24bp total"
echo "  Extracting rightmost 16bp (positions 9-24)"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "DEBUG: Extracting rightmost 16bp from R2..."
    
    # =============================================================================
    # BARCODE EXTRACTION ALGORITHM
    # PURPOSE: Extract rightmost 16bp from each R2 read while preserving FASTQ format
    # 
    # OPERATIONS:
    # - Decompress gzipped FASTQ file using zcat
    # - Process each line based on FASTQ structure (4 lines per read):
    #   * Line 1 (NR%4==1): Header line - keep unchanged
    #   * Line 2 (NR%4==2): Sequence - extract rightmost 16bp using substr()
    #   * Line 3 (NR%4==3): Plus line - keep unchanged
    #   * Line 4 (NR%4==0): Quality scores - extract rightmost 16bp to match sequence
    # - Compress output using gzip
    # 
    # ALGORITHM DETAILS:
    # - substr($0, length($0)-15): Extract from position (length-15) to end
    # - This gives rightmost 16 characters (positions length-15 to length)
    # - Applied to both sequence and quality lines to maintain correspondence
    # 
    # INPUT: Gzipped FASTQ with 24bp R2 reads
    # OUTPUT: Gzipped FASTQ with 16bp extracted barcodes
    # =============================================================================
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) {print $0;}
            else if(NR%4==2) {print substr($0, length($0)-15);}
            else if(NR%4==3) {print $0;}
            else if(NR%4==0) {print substr($0, length($0)-15);}
        }' | gzip > "$EXTRACTED_BC_FILE"
        
    echo "DEBUG: Barcode extraction completed"
else
    echo "DEBUG: Using existing extracted barcodes: $EXTRACTED_BC_FILE"
fi

echo "Extracted barcodes saved to: $EXTRACTED_BC_FILE"

# =============================================================================
# BARCODE EXTRACTION VALIDATION
# PURPOSE: Verify successful barcode extraction and validate output format
# 
# OPERATIONS:
# - Count total number of barcodes (FASTQ lines / 4)
# - Extract sample barcode sequence for length verification
# - Validate that extracted barcodes are exactly 16bp
# 
# VALIDATION CHECKS:
# - Barcode count: Ensures extraction process completed
# - Sample sequence: Provides example of extracted barcode
# - Length validation: Critical check to ensure correct substring extraction
# 
# ERROR HANDLING:
# - Exits with error code 1 if barcode length is not 16bp
# - Prevents downstream analysis with incorrectly extracted barcodes
# 
# OUTPUT EXPECTATIONS:
# - BC_COUNT: Number of successfully extracted barcodes
# - SAMPLE_SEQ: Example 16bp barcode sequence
# - SEQ_LENGTH: Should always be 16 for successful extraction
# =============================================================================
BC_COUNT=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')
SAMPLE_SEQ=$(zcat "$EXTRACTED_BC_FILE" | awk 'NR==2' | head -1)
SEQ_LENGTH=${#SAMPLE_SEQ}

echo "Validation:"
echo "  Total barcodes: $BC_COUNT"
echo "  Sample barcode: $SAMPLE_SEQ"
echo "  Barcode length: $SEQ_LENGTH bp"

if [[ $SEQ_LENGTH -ne 16 ]]; then
    echo "ERROR: Expected 16bp barcodes, got ${SEQ_LENGTH}bp"
    exit 1
fi

# =============================================================================
# COMPLETION SUMMARY
# PURPOSE: Log successful completion of barcode extraction step
# 
# OUTPUT: Timestamp and confirmation of successful processing
# =============================================================================
echo "========================================="
echo "Step 1 complete for $SAMPLE"
echo "End time: $(date)"
echo "=========================================="