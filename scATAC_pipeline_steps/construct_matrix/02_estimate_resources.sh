#!/bin/bash
#SBATCH --job-name=peak_matrix_step2
#SBATCH --output=logs/02_estimate_%a.out
#SBATCH --error=logs/02_estimate_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=30:00
#SBATCH --partition=workq

# ==========================================
# scATAC-seq Pipeline Step 2: Resource Estimation and Disk Space Analysis
# ==========================================
# Purpose: Analyzes input data to estimate computational and storage requirements
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# 
# This script performs comprehensive resource estimation including:
# - Dataset size analysis (peak counts, read counts, overlap estimates)
# - Memory requirements for sorting and bedtools operations
# - Disk space requirements for temporary and output files
# - Processing time estimates based on data size
# - System resource availability checks
# 
# Input: Peaks file, reads file, genome sizes file (via configuration)
# Output: Resource estimation report with recommendations
# Dependencies: 00_config_utils.sh, system utilities (free, df, bc)
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

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps/construct_matrix"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# ==========================================
# STEP 2 FUNCTIONS
# ==========================================

# ==========================================
# DATASET SIZE ESTIMATION FUNCTION
# ==========================================
# Purpose: Analyzes input files to estimate dataset size and complexity
# Implementation: Counts peaks, measures file sizes, estimates read counts and overlaps
# Algorithm:
#   1. Count peaks from peaks file (wc -l)
#   2. Get file sizes using stat command
#   3. Estimate read count from file size (assuming ~50 bytes per compressed line)
#   4. Estimate overlap count using peak-read ratio heuristic
# Input: PEAKS_FILE, READS_FILE, GENOME_SIZES (global variables)
# Output: Exported estimation variables for downstream use
# Edge cases: Handles compressed files, provides conservative estimates
# Returns: 0 on success
estimate_dataset_size() {
    log_step "DATASET_ESTIMATION" "START"
    logger "Estimating dataset size for sample: $SAMPLE"
    
    # Get file sizes and counts
    local peak_count=$(wc -l < "$PEAKS_FILE")
    local reads_size=$(stat -c%s "$READS_FILE")
    local genome_size=$(stat -c%s "$GENOME_SIZES")
    
    # Estimate compressed reads count (rough approximation)
    local estimated_reads=$((reads_size / 50))  # Assuming ~50 bytes per compressed read line
    
    # Estimate overlap records (conservative estimate)
    local estimated_overlaps=$((peak_count * estimated_reads / 1000000))
    
    logger "Dataset estimation:"
    logger "  Peak count: $(printf "%'d" $peak_count)"
    logger "  Reads file size: $(numfmt --to=iec-i --suffix=B $reads_size)"
    logger "  Estimated reads: $(printf "%'d" $estimated_reads)"
    logger "  Estimated overlap records: $(printf "%'d" $estimated_overlaps)"
    
    # Export estimates for use in other steps
    export ESTIMATED_PEAK_COUNT=$peak_count
    export ESTIMATED_READS_SIZE=$reads_size
    export ESTIMATED_READS_COUNT=$estimated_reads
    export ESTIMATED_OVERLAPS=$estimated_overlaps
    
    log_step "DATASET_ESTIMATION" "COMPLETE"
    return 0
}

# ==========================================
# MEMORY REQUIREMENTS ESTIMATION FUNCTION
# ==========================================
# Purpose: Calculates memory requirements for each pipeline step
# Implementation: Estimates memory based on file sizes and data complexity
# Algorithm:
#   1. Base memory requirements for different operations
#   2. Dynamic scaling based on estimated overlap count
#   3. System memory availability check
#   4. Conservative estimates with safety margins
# Input: ESTIMATED_OVERLAPS (from dataset estimation)
# Output: Exported memory recommendation variables
# Edge cases: Provides warnings for insufficient system memory
# Returns: 0 on success
estimate_memory_requirements() {
    log_step "MEMORY_ESTIMATION" "START"
    logger "Estimating memory requirements..."
    
    # Base memory requirements
    local base_memory_gb=8
    local sort_memory_gb=32
    local bedtools_memory_gb=16
    
    # Adjust based on dataset size
    if [[ $ESTIMATED_OVERLAPS -gt 100000000 ]]; then
        sort_memory_gb=64
        bedtools_memory_gb=32
        logger "Large dataset detected, increasing memory estimates"
    elif [[ $ESTIMATED_OVERLAPS -gt 50000000 ]]; then
        sort_memory_gb=48
        bedtools_memory_gb=24
        logger "Medium dataset detected, moderate memory estimates"
    fi
    
    # Get current system memory
    local total_memory_gb=0
    local available_memory_gb=0
    
    if command -v free &> /dev/null; then
        total_memory_gb=$(free -g | awk '/^Mem:/ {print $2}')
        available_memory_gb=$(free -g | awk '/^Mem:/ {print $7}')
    fi
    
    logger "Memory estimation:"
    logger "  System total memory: ${total_memory_gb}GB"
    logger "  System available memory: ${available_memory_gb}GB"
    logger "  Recommended sort memory: ${sort_memory_gb}GB"
    logger "  Recommended bedtools memory: ${bedtools_memory_gb}GB"
    
    # Check if system has sufficient memory
    if [[ $available_memory_gb -lt $sort_memory_gb ]]; then
        logger "WARNING: Available memory (${available_memory_gb}GB) is less than recommended sort memory (${sort_memory_gb}GB)"
        logger "Consider using a node with more memory or reducing sort memory"
    fi
    
    # Export memory estimates
    export RECOMMENDED_SORT_MEMORY="${sort_memory_gb}G"
    export RECOMMENDED_BEDTOOLS_MEMORY="${bedtools_memory_gb}G"
    
    log_step "MEMORY_ESTIMATION" "COMPLETE"
    return 0
}

# ==========================================
# DISK REQUIREMENTS ESTIMATION FUNCTION
# ==========================================
# Purpose: Calculates disk space requirements for all pipeline stages
# Implementation: Estimates space for input, temporary, and output files with safety margins
# Algorithm:
#   1. Calculate input data size from reads file
#   2. Estimate temporary space (3x input for sorting and processing)
#   3. Estimate output space (1x input for final matrices)
#   4. Add safety buffer and check against available space
#   5. Generate warnings for insufficient space
# Input: ESTIMATED_READS_SIZE, TEMP_DIR, OUTPUT_DIR (global variables)
# Output: Exported disk space variables and availability warnings
# Edge cases: Handles missing directories, provides conservative estimates
# Returns: 0 on success
estimate_disk_requirements() {
    log_step "DISK_ESTIMATION" "START"
    logger "Estimating disk space requirements..."
    
    # Base disk requirements (in GB)
    local input_size_gb=$((ESTIMATED_READS_SIZE / 1024 / 1024 / 1024))
    local temp_multiplier=3  # Temporary files can be 3x input size
    local output_multiplier=1  # Output files roughly same size as input
    
    local temp_disk_gb=$((input_size_gb * temp_multiplier))
    local output_disk_gb=$((input_size_gb * output_multiplier))
    local total_required_gb=$((temp_disk_gb + output_disk_gb + 50))  # 50GB buffer
    
    # Check current disk space
    local temp_available_gb=0
    local output_available_gb=0
    
    if [[ -d "$TEMP_DIR" ]]; then
        local temp_available_kb=$(df -k "$TEMP_DIR" | awk 'NR==2 {print $4}')
        temp_available_gb=$((temp_available_kb / 1024 / 1024))
    fi
    
    if [[ -d "$OUTPUT_DIR" ]]; then
        local output_available_kb=$(df -k "$OUTPUT_DIR" | awk 'NR==2 {print $4}')
        output_available_gb=$((output_available_kb / 1024 / 1024))
    fi
    
    logger "Disk space estimation:"
    logger "  Input data size: ${input_size_gb}GB"
    logger "  Estimated temp space needed: ${temp_disk_gb}GB"
    logger "  Estimated output space needed: ${output_disk_gb}GB"
    logger "  Total space recommended: ${total_required_gb}GB"
    logger "  Temp directory available: ${temp_available_gb}GB"
    logger "  Output directory available: ${output_available_gb}GB"
    
    # Check if sufficient disk space is available
    local warnings=()
    
    if [[ $temp_available_gb -lt $temp_disk_gb ]]; then
        warnings+=("Insufficient temp space: ${temp_available_gb}GB available, ${temp_disk_gb}GB needed")
    fi
    
    if [[ $output_available_gb -lt $output_disk_gb ]]; then
        warnings+=("Insufficient output space: ${output_available_gb}GB available, ${output_disk_gb}GB needed")
    fi
    
    if [[ ${#warnings[@]} -gt 0 ]]; then
        logger "DISK SPACE WARNINGS:"
        for warning in "${warnings[@]}"; do
            logger "  WARNING: $warning"
        done
    else
        logger "Disk space check: PASSED"
    fi
    
    # Export disk estimates
    export ESTIMATED_TEMP_DISK_GB=$temp_disk_gb
    export ESTIMATED_OUTPUT_DISK_GB=$output_disk_gb
    export ESTIMATED_TOTAL_DISK_GB=$total_required_gb
    
    log_step "DISK_ESTIMATION" "COMPLETE"
    return 0
}

# ==========================================
# PROCESSING TIME ESTIMATION FUNCTION
# ==========================================
# Purpose: Estimates processing time for each pipeline step
# Implementation: Uses base time estimates with scaling factors based on data size
# Algorithm:
#   1. Define base processing times for bedtools, overlap, and matrix steps
#   2. Calculate size factor based on estimated overlap count (simple thresholds)
#   3. Scale individual step times by size factor
#   4. Sum total processing time and convert to hours using bc
# Input: ESTIMATED_OVERLAPS (from dataset estimation)
# Output: Exported time estimation variables for job scheduling
# Edge cases: Handles bc command availability, provides fallback calculations
# Returns: 0 on success
estimate_processing_time() {
    log_step "TIME_ESTIMATION" "START"
    logger "Estimating processing time..."
    
    # Base time estimates (in minutes)
    local bedtools_base_minutes=30
    local overlap_processing_minutes=20
    local matrix_creation_minutes=15
    
    # Adjust based on dataset size
    local size_factor=1
    if [[ $ESTIMATED_OVERLAPS -gt 100000000 ]]; then
        size_factor=4
    elif [[ $ESTIMATED_OVERLAPS -gt 50000000 ]]; then
        size_factor=2
    fi
    
    local estimated_bedtools_minutes=$((bedtools_base_minutes * size_factor))
    local estimated_overlap_minutes=$((overlap_processing_minutes * size_factor))
    local estimated_matrix_minutes=$((matrix_creation_minutes * size_factor))
    local total_minutes=$((estimated_bedtools_minutes + estimated_overlap_minutes + estimated_matrix_minutes))
    
    logger "Processing time estimation:"
    logger "  Bedtools intersect: ${estimated_bedtools_minutes} minutes"
    logger "  Overlap processing: ${estimated_overlap_minutes} minutes"
    logger "  Matrix creation: ${estimated_matrix_minutes} minutes"
    logger "  Total estimated time: ${total_minutes} minutes ($(echo "scale=1; $total_minutes/60" | bc -l 2>/dev/null || echo "N/A") hours)"
    
    # Export time estimates
    export ESTIMATED_BEDTOOLS_MINUTES=$estimated_bedtools_minutes
    export ESTIMATED_OVERLAP_MINUTES=$estimated_overlap_minutes
    export ESTIMATED_MATRIX_MINUTES=$estimated_matrix_minutes
    export ESTIMATED_TOTAL_MINUTES=$total_minutes
    
    log_step "TIME_ESTIMATION" "COMPLETE"
    return 0
}

# ==========================================
# RESOURCE REPORT GENERATION FUNCTION
# ==========================================
# Purpose: Creates comprehensive resource estimation report
# Implementation: Generates formatted text report with all estimation results
# Algorithm:
#   1. Create report file with timestamp and job information
#   2. Include dataset size analysis results
#   3. Add memory and disk space requirements
#   4. Include processing time estimates
#   5. Provide recommendations for job submission
# Input: All exported estimation variables from previous functions
# Output: Formatted text report file in output directory
# Edge cases: Handles missing system information gracefully
# Returns: Implicit (void function)
generate_resource_report() {
    local report_file="$OUTPUT_DIR/step2_resource_estimation_${SAMPLE}.txt"
    
    logger "Generating resource estimation report: $report_file"
    
    cat > "$report_file" << EOF
STEP 2 RESOURCE ESTIMATION REPORT
==================================
Sample: $SAMPLE
Estimation Time: $(date)
SLURM Job ID: ${SLURM_JOB_ID:-"N/A"}

DATASET SIZE:
- Peak count: $(printf "%'d" $ESTIMATED_PEAK_COUNT)
- Reads file size: $(numfmt --to=iec-i --suffix=B $ESTIMATED_READS_SIZE)
- Estimated reads: $(printf "%'d" $ESTIMATED_READS_COUNT)
- Estimated overlaps: $(printf "%'d" $ESTIMATED_OVERLAPS)

MEMORY REQUIREMENTS:
- Recommended sort memory: $RECOMMENDED_SORT_MEMORY
- Recommended bedtools memory: $RECOMMENDED_BEDTOOLS_MEMORY
- System total memory: $(free -g 2>/dev/null | awk '/^Mem:/ {print $2"GB"}' || echo "N/A")
- System available memory: $(free -g 2>/dev/null | awk '/^Mem:/ {print $7"GB"}' || echo "N/A")

DISK SPACE REQUIREMENTS:
- Estimated temp space needed: ${ESTIMATED_TEMP_DISK_GB}GB
- Estimated output space needed: ${ESTIMATED_OUTPUT_DISK_GB}GB
- Total space recommended: ${ESTIMATED_TOTAL_DISK_GB}GB
- Temp directory available: $(df -h "$TEMP_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo "N/A")
- Output directory available: $(df -h "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2 {print $4}' || echo "N/A")

PROCESSING TIME ESTIMATES:
- Bedtools intersect: ${ESTIMATED_BEDTOOLS_MINUTES} minutes
- Overlap processing: ${ESTIMATED_OVERLAP_MINUTES} minutes
- Matrix creation: ${ESTIMATED_MATRIX_MINUTES} minutes
- Total estimated time: ${ESTIMATED_TOTAL_MINUTES} minutes

RECOMMENDATIONS:
- Use SLURM job with at least $(echo "scale=0; $ESTIMATED_TOTAL_MINUTES * 1.5 / 60" | bc -l 2>/dev/null || echo "8") hours time limit
- Request at least $(echo "$RECOMMENDED_SORT_MEMORY" | sed 's/G//') GB memory
- Ensure ${ESTIMATED_TOTAL_DISK_GB}GB disk space is available
- Monitor temp directory usage during processing
EOF

    logger "Step 2 resource estimation report saved: $report_file"
}

# ==========================================
# MAIN EXECUTION
# ==========================================
# Purpose: Orchestrates the complete resource estimation workflow
# Implementation: Calls all estimation functions in sequence and generates report
# Algorithm:
#   1. Initialize logging and timing for the step
#   2. Setup temporary directory if needed
#   3. Run dataset size estimation with progress tracking
#   4. Calculate memory requirements
#   5. Estimate disk space needs
#   6. Calculate processing time
#   7. Generate comprehensive report
#   8. Log completion status and timing
# Input: Sample configuration from global variables
# Output: Resource estimation report file and exported variables
# Edge cases: Handles failures in individual estimation steps
# Returns: 0 on success, exits on critical failures

main() {
    local start_time=$(date +%s)
    
    # Initialize logging at the very beginning
    # This ensures LOG_FILE is defined before any logger calls
    init_logging

    echo "========================================="
    echo "Step 2: Resource Estimation"
    echo "Sample: $SAMPLE"
    echo "Start time: $(date)"
    echo "SLURM Job ID: ${SLURM_JOB_ID:-"N/A"}"
    echo "========================================="
    
    # Setup temporary directory if not already done
    if [[ ! -d "$TEMP_DIR" ]]; then
        setup_temp_dir
    fi
    
    # Estimate dataset size
    log_progress "estimation" 0
    estimate_dataset_size || {
        log_step "STEP2" "FAILED"
        exit 1
    }
    log_progress "estimation" 25
    
    # Estimate memory requirements
    estimate_memory_requirements || {
        log_step "STEP2" "FAILED"
        exit 1
    }
    log_progress "estimation" 50
    
    # Estimate disk requirements
    estimate_disk_requirements || {
        log_step "STEP2" "FAILED"
        exit 1
    }
    log_progress "estimation" 75
    
    # Estimate processing time
    estimate_processing_time || {
        log_step "STEP2" "FAILED"
        exit 1
    }
    log_progress "estimation" 90
    
    # Generate resource estimation report
    generate_resource_report
    log_progress "estimation" 100
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_step "STEP2" "COMPLETE"
    
    echo "========================================="
    echo "Step 2 completed successfully for $SAMPLE"
    echo "Processing time: $(date -d@$duration -u +%H:%M:%S)"
    echo "Resource report: $OUTPUT_DIR/step2_resource_estimation_${SAMPLE}.txt"
    echo "Log file: $LOG_FILE"
    echo "Ready for Step 3: Bedtools Intersect"
    echo "End time: $(date)"
    echo "========================================="
}

# ==========================================
# SCRIPT EXECUTION ENTRY POINT
# ==========================================
# Purpose: Entry point for script execution with error handling
# Implementation: Calls main function with all command line arguments
# Input: Command line arguments (passed through to main)
# Output: Script execution results and exit codes
# Edge cases: Ensures proper error propagation from main function

# Execute main function
main "$@"