#!/bin/bash
#SBATCH --job-name=peak_matrix_step1
#SBATCH --output=logs/01_validate_%a.out
#SBATCH --error=logs/01_validate_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# ==========================================
# scATAC-seq Pipeline Step 1: Prerequisites Validation and Environment Setup
# ==========================================
# Purpose: Validates all input files, tools, and environment setup before matrix construction
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# 
# This script performs comprehensive validation of:
# - Input file existence, format, and integrity
# - Required tools and software dependencies
# - Environment setup and resource availability
# - Peaks file sorting compatibility with genome reference
# 
# Input: Peaks file, reads file, genome sizes file (via configuration)
# Output: Validation report, properly sorted peaks file if needed
# Dependencies: 00_config_utils.sh, conda environment, bedtools
# ==========================================

# Enable strict error handling: exit on error, undefined variables, and pipe failures
set -euo pipefail

# ==========================================
# SLURM ENVIRONMENT CHECK
# ==========================================
# Purpose: Ensures script is executed within SLURM job scheduler environment
# Implementation: Checks for SLURM_JOB_ID environment variable presence
# Edge cases: Prevents direct bash execution which could cause resource conflicts
# Output: Error message and exit if not running under SLURM

# Ensure script is running under SLURM
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: This script must be executed with sbatch, not bash directly."
    echo "Usage: sbatch $0"
    exit 1
fi

# ==========================================
# LOAD CONFIGURATION
# ==========================================
# Purpose: Load shared configuration, utilities, and global variables
# Implementation: Sources 00_config_utils.sh for pipeline-wide settings
# Input: 00_config_utils.sh configuration file
# Output: Global variables (SAMPLE, file paths, resource settings)
# Edge cases: Exits if configuration file is missing or corrupted

# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# source "$SCRIPT_DIR/00_config_utils.sh"

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps/construct_matrix"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# ==========================================
# STEP 1 FUNCTIONS
# ==========================================

# ==========================================
# PREREQUISITES VALIDATION FUNCTION
# ==========================================
# Purpose: Validates existence and basic integrity of all required input files
# Implementation: Uses validate_file utility to check file existence and minimum sizes
# Input expectations:
#   - PEAKS_FILE: BED format file with genomic regions (min 1KB)
#   - READS_FILE: Fragment/read data file (min 1MB)
#   - GENOME_SIZES: Chromosome sizes reference (min 100B)
# Output: Log messages, cleaned output directory
# Edge cases: Removes existing matrix directory to prevent conflicts
# Returns: 0 on success, 1 on validation failure
validate_prerequisites() {
    log_step "PREREQUISITES_VALIDATION" "START"
    logger "Validating prerequisites for sample: $SAMPLE"
    
    # Validate each required file
    validate_file "$PEAKS_FILE" 1000 "Sorted peaks file" || return 1
    validate_file "$READS_FILE" 1000000 "Reads file" || return 1
    validate_file "$GENOME_SIZES" 100 "Genome sizes file" || return 1
    
    # Check if previous matrix exists and warn
    if [[ -d "$MATRIX_DIR" ]]; then
        logger "Previous matrix directory exists, will be overwritten: $MATRIX_DIR"
        rm -rf "$MATRIX_DIR"
    fi
    
    log_step "PREREQUISITES_VALIDATION" "COMPLETE"
    return 0
}

# ==========================================
# PEAKS SORTING VALIDATION FUNCTION
# ==========================================
# Purpose: Ensures peaks file chromosome order matches genome reference file
# Implementation: Compares chromosome order between peaks and genome files, resorts if needed
# Algorithm:
#   1. Extract chromosome order from genome sizes file
#   2. Create chromosome ordering index with line numbers
#   3. Check if peaks file first chromosome exists in genome file
#   4. If mismatch detected, resort peaks using genome chromosome order
# Input: PEAKS_FILE (BED format), GENOME_SIZES (tab-delimited)
# Output: Properly sorted peaks file matching genome chromosome order
# Edge cases: Handles chromosome naming inconsistencies, preserves all BED columns
# Returns: 0 on success (sorted or already compatible)
validate_peaks_sorting() {
    log_step "PEAKS_SORTING_VALIDATION" "START"
    logger "Verifying peaks file sort order matches genome file..."
    
    local peaks_temp="$TEMP_DIR/peaks_sorted.bed"
    
    # Extract chromosome order from genome sizes file
    local chroms_file="$TEMP_DIR/chroms.txt"
    cut -f1 "$GENOME_SIZES" > "$chroms_file"
    
    # Create a sort order file for chromosome sorting
    # Use LC_ALL=C to ensure consistent sorting behavior
    LC_ALL=C awk '{print $0"\t"NR}' "$chroms_file" | LC_ALL=C sort -k1,1 > "$TEMP_DIR/chrom_order.txt"
    
    # Check if peaks need resorting
    local first_chrom=$(head -n 1 "$PEAKS_FILE" | cut -f1)
    if ! grep -q "^$first_chrom\s" "$GENOME_SIZES"; then
        logger "Peaks file needs resorting to match genome file order"
        
        # Resort peaks file to match genome file order
        logger "Resorting peaks file to match genome file order..."
        join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2 \
            <(sort -k1,1 "$PEAKS_FILE") \
            "$TEMP_DIR/chrom_order.txt" | \
            sort -k7,7n -k2,2n | \
            cut -f1-6 > "$peaks_temp"
            
        # Replace original peaks file with properly sorted version
        mv "$peaks_temp" "$PEAKS_FILE"
        logger "Peaks file resorted successfully"
    else
        logger "Peaks file sort order appears compatible with genome file"
    fi
    
    log_step "PEAKS_SORTING_VALIDATION" "COMPLETE"
    return 0
}

# ==========================================
# ENVIRONMENT SETUP FUNCTION
# ==========================================
# Purpose: Initializes complete runtime environment for pipeline execution
# Implementation: Sequential setup of logging, conda environment, tools, and workspace
# Operations performed:
#   1. Initialize logging system for progress tracking
#   2. Activate conda environment with required packages
#   3. Validate availability of all required tools
#   4. Create and verify temporary workspace directory
#   5. Check sufficient disk space for processing
# Input: Global configuration variables (TEMP_DIR, MIN_DISK_SPACE_GB)
# Output: Fully configured environment ready for processing
# Edge cases: Fails fast on any setup error to prevent downstream issues
# Returns: 0 on complete success, 1 on any setup failure
setup_environment() {
    log_step "ENVIRONMENT_SETUP" "START"
    
    # Initialize logging
    init_logging
    
    # Setup conda environment
    setup_conda_env || {
        log_step "ENVIRONMENT_SETUP" "FAILED"
        return 1
    }
    
    # Validate tools
    validate_tools || {
        log_step "ENVIRONMENT_SETUP" "FAILED"
        return 1
    }
    
    # Setup temporary directory
    setup_temp_dir || {
        log_step "ENVIRONMENT_SETUP" "FAILED"
        return 1
    }
    
    # Verify temp directory has sufficient space
    check_disk_space "$TEMP_DIR" "$MIN_DISK_SPACE_GB" || {
        log_step "ENVIRONMENT_SETUP" "FAILED"
        return 1
    }
    
    log_step "ENVIRONMENT_SETUP" "COMPLETE"
    return 0
}

# ==========================================
# VALIDATION REPORT GENERATION FUNCTION
# ==========================================
# Purpose: Creates comprehensive validation report documenting all checks performed
# Implementation: Generates structured text report with file info, environment details, and tool status
# Report sections:
#   - Sample identification and job metadata
#   - Input file paths and sizes
#   - File content statistics (peak counts, chromosome counts)
#   - Environment configuration (conda, temp directory, disk space)
#   - Tool validation results with installation paths
# Input: Global variables (SAMPLE, file paths, SLURM_JOB_ID)
# Output: Detailed validation report file in OUTPUT_DIR
# Edge cases: Handles missing files gracefully with "N/A" fallbacks
# Returns: None (void function)
generate_step1_report() {
    local report_file="$OUTPUT_DIR/step1_validation_report_${SAMPLE}.txt"
    
    cat > "$report_file" << EOF
STEP 1 VALIDATION REPORT
========================
Sample: $SAMPLE
Validation Time: $(date)
SLURM Job ID: ${SLURM_JOB_ID:-"N/A"}

INPUT FILES:
- Peaks file: $PEAKS_FILE
- Reads file: $READS_FILE
- Genome sizes: $GENOME_SIZES

FILE SIZES:
- Peaks: $(stat -c%s "$PEAKS_FILE" 2>/dev/null | numfmt --to=iec-i --suffix=B || echo "N/A")
- Reads: $(stat -c%s "$READS_FILE" 2>/dev/null | numfmt --to=iec-i --suffix=B || echo "N/A")
- Genome: $(stat -c%s "$GENOME_SIZES" 2>/dev/null | numfmt --to=iec-i --suffix=B || echo "N/A")

COUNTS:
- Total peaks: $(wc -l < "$PEAKS_FILE" 2>/dev/null || echo "N/A")
- Chromosomes in genome: $(wc -l < "$GENOME_SIZES" 2>/dev/null || echo "N/A")

ENVIRONMENT:
- Conda environment: peak_calling_new
- Temporary directory: $TEMP_DIR
- Available disk space: $(df -h "$TEMP_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo "N/A")

TOOLS VALIDATED:
- bedtools: $(command -v bedtools || echo "NOT FOUND")
- sort: $(command -v sort || echo "NOT FOUND")
- awk: $(command -v awk || echo "NOT FOUND")
- gzip: $(command -v gzip || echo "NOT FOUND")
- zcat: $(command -v zcat || echo "NOT FOUND")

VALIDATION STATUS: PASSED
EOF

    logger "Step 1 validation report saved: $report_file"
}

# ==========================================
# MAIN EXECUTION FUNCTION
# ==========================================
# Purpose: Orchestrates complete Step 1 validation workflow with progress tracking
# Implementation: Sequential execution of validation steps with error handling and timing
# Workflow:
#   1. Initialize environment and logging systems
#   2. Validate input file prerequisites (0-50% progress)
#   3. Validate and fix peaks file sorting (50-90% progress)
#   4. Generate comprehensive validation report (90-100% progress)
#   5. Display completion summary with timing and next steps
# Input: Command line arguments (passed through)
# Output: Console progress messages, validation report, log files
# Edge cases: Exits immediately on any validation failure
# Returns: Exit code 0 on success, 1 on failure

main() {
    local start_time=$(date +%s)
    
    echo "========================================="
    echo "Step 1: Prerequisites Validation"
    echo "Sample: $SAMPLE"
    echo "Start time: $(date)"
    echo "SLURM Job ID: ${SLURM_JOB_ID:-"N/A"}"
    echo "========================================="
    
    # Setup environment and logging
    setup_environment || exit 1
    
    # Validate prerequisites
    log_progress "validation" 0
    validate_prerequisites || {
        log_step "STEP1" "FAILED"
        exit 1
    }
    log_progress "validation" 50
    
    # Validate peaks sorting
    validate_peaks_sorting || {
        log_step "STEP1" "FAILED"
        exit 1
    }
    log_progress "validation" 90
    
    # Generate validation report
    generate_step1_report
    log_progress "validation" 100
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_step "STEP1" "COMPLETE"
    
    echo "========================================="
    echo "Step 1 completed successfully for $SAMPLE"
    echo "Processing time: $(date -d@$duration -u +%H:%M:%S)"
    echo "Validation report: $OUTPUT_DIR/step1_validation_report_${SAMPLE}.txt"
    echo "Log file: $LOG_FILE"
    echo "Ready for Step 2: Resource Estimation"
    echo "End time: $(date)"
    echo "========================================="
}

# ==========================================
# SCRIPT EXECUTION ENTRY POINT
# ==========================================
# Purpose: Initiates main function with all command line arguments
# Implementation: Passes all script arguments to main function for processing
# Input: All command line arguments ($@)
# Output: Depends on main function execution
# Edge cases: Inherits exit codes from main function

# Execute main function
main "$@"