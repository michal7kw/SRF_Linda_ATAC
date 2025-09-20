#!/bin/bash
# =============================================================================
# SCRIPT: 03_validate_counts.sh
# PURPOSE: Validate read count consistency across FASTQ files in scATAC-seq pipeline
# 
# DESCRIPTION:
# This script ensures that R1 (forward genomic), R3 (reverse genomic), and 
# extracted R2 barcode files contain the same number of reads. Mismatched counts
# indicate data corruption, incomplete files, or extraction errors that would
# cause downstream alignment failures.
# 
# KEY OPERATIONS:
# 1. Count reads in R1, R3, and extracted barcode files
# 2. Validate count consistency across all files
# 3. Attempt automatic repair if barcode extraction failed
# 4. Generate validation report for pipeline tracking
# 
# INPUT REQUIREMENTS:
# - Original R1, R2, R3 FASTQ files from sequencing
# - Extracted barcode file from step 1 (01_extract_barcodes.sh)
# 
# OUTPUT:
# - Validation status report
# - Corrected barcode file (if repair was needed)
# - Error exit if validation fails
# 
# EDGE CASES:
# - Handles barcode extraction failures with automatic re-extraction
# - Detects and reports corrupted or incomplete FASTQ files
# - Validates original R2 file integrity before attempting repair
# =============================================================================
#SBATCH --job-name=validate_counts
#SBATCH --output=logs/03_validate_counts_%a.out
#SBATCH --error=logs/03_validate_counts_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# =============================================================================
# ENVIRONMENT SETUP
# PURPOSE: Initialize conda environment and error handling
# 
# OPERATIONS:
# - Activates alignment_two conda environment with required tools
# - Enables strict error handling (exit on error, undefined variables, pipe failures)
# 
# REQUIRED TOOLS:
# - Standard Unix utilities (zcat, wc, awk, gzip)
# - File system operations for validation and repair
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# SAMPLE AND PATH CONFIGURATION
# PURPOSE: Define sample identifiers and directory structure
# 
# SAMPLE ARRAY:
# - Uses SLURM array indexing for parallel processing
# - Each array job processes one sample independently
# 
# DIRECTORY STRUCTURE:
# - DATA_DIR: Original FASTQ files from sequencing
# - PROCESSED_DATA_DIR: Intermediate files (extracted barcodes)
# - OUTPUT_DIR: Final pipeline outputs and QC reports
# 
# FILE NAMING CONVENTION:
# - R1: Forward genomic reads (${SAMPLE}_R1_001.fastq.gz)
# - R2: Barcode reads (${SAMPLE}_R2_001.fastq.gz)
# - R3: Reverse genomic reads (${SAMPLE}_R3_001.fastq.gz)
# - Extracted: Processed barcodes (${SAMPLE}_R2_rightmost16bp.fastq.gz)
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 3: Validating read counts for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# INPUT FILE VALIDATION
# PURPOSE: Verify that prerequisite files exist before validation
# 
# DEPENDENCIES:
# - Extracted barcode file from step 1 (01_extract_barcodes.sh)
# - Original R1, R2, R3 FASTQ files from sequencing
# 
# ERROR HANDLING:
# - Exits with clear error message if barcode file missing
# - Provides guidance on required prerequisite steps
# =============================================================================
EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run step 1 (01_extract_barcodes.sh) first"
    exit 1
fi

# =============================================================================
# READ COUNT ANALYSIS
# PURPOSE: Count reads in all FASTQ files for consistency validation
# 
# ALGORITHM:
# - Uses zcat to decompress gzipped FASTQ files
# - Counts total lines and divides by 4 (FASTQ format: 4 lines per read)
# - Processes each file independently for parallel efficiency
# 
# FASTQ FORMAT:
# - Line 1: Sequence identifier (@...)
# - Line 2: Raw sequence
# - Line 3: Quality identifier (+...)
# - Line 4: Quality scores
# 
# EXPECTED BEHAVIOR:
# - All three files should contain identical read counts
# - Mismatches indicate data corruption or processing errors
# =============================================================================
echo "DEBUG: Counting reads in each file..."
R1_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" | wc -l | awk '{print $1/4}')
R3_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" | wc -l | awk '{print $1/4}') 
BC_COUNT=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')

echo "Read counts:"
echo "  R1 (genomic, forward): $R1_COUNT"
echo "  R3 (genomic, reverse): $R3_COUNT" 
echo "  R2 (extracted barcodes): $BC_COUNT"

# =============================================================================
# PRIMARY VALIDATION: R1 vs R3 COUNT CONSISTENCY
# PURPOSE: Verify that forward and reverse genomic reads have matching counts
# 
# RATIONALE:
# - R1 and R3 are paired-end genomic reads from the same sequencing run
# - Count mismatch indicates fundamental data corruption or incomplete files
# - This is a critical error that cannot be automatically corrected
# 
# ERROR HANDLING:
# - Immediate exit if counts don't match
# - Requires manual investigation of original FASTQ files
# =============================================================================
if [[ "$R1_COUNT" -ne "$R3_COUNT" ]]; then
    echo "ERROR: R1 and R3 read counts don't match ($R1_COUNT vs $R3_COUNT)"
    exit 1
fi

# =============================================================================
# BARCODE COUNT VALIDATION AND AUTOMATIC REPAIR
# PURPOSE: Validate barcode file consistency and attempt repair if needed
# 
# VALIDATION LOGIC:
# - Compares extracted barcode count with R1 genomic read count
# - Mismatches suggest barcode extraction errors (not data corruption)
# 
# AUTOMATIC REPAIR STRATEGY:
# 1. Verify original R2 file integrity by comparing with R1 count
# 2. Re-extract barcodes using same algorithm as step 1
# 3. Validate that repair was successful
# 
# ERROR CONDITIONS:
# - Original R2 file has different count than R1 (data corruption)
# - Re-extraction fails to produce matching counts (algorithm error)
# 
# REPAIR ALGORITHM:
# - Extracts rightmost 16bp from R2 reads (positions 9-24)
# - Maintains FASTQ format with original headers and quality scores
# - Uses awk for efficient line-by-line processing
# =============================================================================
if [[ "$R1_COUNT" -ne "$BC_COUNT" ]]; then
    echo "ERROR: R1 and barcode counts don't match ($R1_COUNT vs $BC_COUNT)"
    echo "This suggests the R2 barcode file was truncated or has different structure"
    echo "Attempting to fix by re-extracting barcodes..."
    
    # Re-extract barcodes ensuring same read count as R1
    echo "DEBUG: Re-extracting barcodes with read count validation..."
    ORIG_R2_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | wc -l | awk '{print $1/4}')
    echo "DEBUG: Original R2 read count: $ORIG_R2_COUNT"
    
    if [[ "$ORIG_R2_COUNT" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Original R2 file has different read count than R1 ($ORIG_R2_COUNT vs $R1_COUNT)"
        echo "This indicates corrupted or incomplete FASTQ files"
        exit 1
    fi
    
    # =============================================================================
    # BARCODE RE-EXTRACTION PROCESS
    # PURPOSE: Regenerate barcode file with correct read count
    # 
    # EXTRACTION ALGORITHM:
    # - Line 1 (NR%4==1): Copy sequence header unchanged
    # - Line 2 (NR%4==2): Extract substring from position 9, length 16
    # - Line 3 (NR%4==3): Copy quality header unchanged  
    # - Line 4 (NR%4==0): Extract corresponding quality scores (positions 9-24)
    # 
    # QUALITY CONTROL:
    # - Removes existing file to prevent partial overwrites
    # - Compresses output to maintain consistency with pipeline
    # - Validates extraction success before proceeding
    # 
    # POSITION RATIONALE:
    # - Position 9-24: Rightmost 16bp of 24bp barcode sequence
    # - Matches extraction logic from step 1 (01_extract_barcodes.sh)
    # =============================================================================
    rm -f "$EXTRACTED_BC_FILE"
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,9,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,9,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
    
    # Recount after extraction
    BC_COUNT_NEW=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')
    echo "DEBUG: New barcode count after re-extraction: $BC_COUNT_NEW"
    
    if [[ "$BC_COUNT_NEW" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Re-extraction failed, counts still don't match ($BC_COUNT_NEW vs $R1_COUNT)"
        exit 1
    fi
    
    echo "DEBUG: Re-extraction successful, all counts now match"
    BC_COUNT=$BC_COUNT_NEW
fi

# =============================================================================
# VALIDATION SUCCESS CONFIRMATION
# PURPOSE: Confirm that all files have consistent read counts
# 
# SUCCESS CRITERIA:
# - R1, R3, and extracted barcode files all contain identical read counts
# - No data corruption or extraction errors detected
# - Pipeline can proceed safely to alignment step
# =============================================================================
echo "✓ All file read counts match: $R1_COUNT reads per file"

# =============================================================================
# VALIDATION RESULTS PERSISTENCE
# PURPOSE: Save validation results for pipeline tracking and debugging
# 
# OPERATIONS:
# - Creates QC directory if it doesn't exist
# - Generates structured validation report
# - Records all read counts and validation status
# - Timestamps validation for audit trail
# 
# REPORT CONTENTS:
# - Sample identifier and read counts for all files
# - Validation status (PASS/FAIL)
# - Timestamp for tracking when validation occurred
# 
# USAGE:
# - Downstream scripts can verify validation was completed
# - Debugging tool for investigating count mismatches
# - Quality control documentation for pipeline runs
# =============================================================================
mkdir -p "$OUTPUT_DIR/qc"
cat > "$OUTPUT_DIR/qc/${SAMPLE}_read_count_validation.txt" << EOF
SAMPLE=$SAMPLE
R1_COUNT=$R1_COUNT
R3_COUNT=$R3_COUNT
BC_COUNT=$BC_COUNT
VALIDATION_STATUS=PASS
VALIDATED_AT=$(date)
EOF

echo "Validation results saved to: $OUTPUT_DIR/qc/${SAMPLE}_read_count_validation.txt"

# =============================================================================
# COMPLETION SUMMARY
# PURPOSE: Log successful validation completion with key metrics
# 
# OUTPUT: Confirmation of validation success and final read count
# =============================================================================
echo "========================================="
echo "Step 3 complete for $SAMPLE"
echo "All read counts validated: $R1_COUNT reads"
echo "End time: $(date)"
echo "=========================================="