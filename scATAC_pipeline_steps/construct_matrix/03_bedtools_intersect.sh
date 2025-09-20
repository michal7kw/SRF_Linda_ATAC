#!/bin/bash
#SBATCH --job-name=03opt_bedtools_intersect
#SBATCH --output=logs/03_bedtools_intersect_opt_%a.out
#SBATCH --error=logs/03_bedtools_intersect_opt_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq

# ==========================================
# scATAC-seq Pipeline Step 3: Optimized Bedtools Intersect Operation
# ==========================================
# Purpose: Performs optimized intersection between peaks and fragment files
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# 
# This script performs high-performance intersection analysis including:
# - Optimized sorting of fragment files with parallel processing
# - Efficient bedtools intersect operation with sorted inputs
# - Memory-optimized temporary file handling
# - Comprehensive validation and error handling
# - Performance monitoring and statistics generation
# 
# Input: Peaks file (sorted), fragments file (BED format)
# Output: Intersection results file with peak-fragment overlaps
# Dependencies: 00_config_utils.sh, bedtools, system utilities (sort, wc)
# Performance: Uses parallel processing and optimized I/O for large datasets
# ==========================================

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
# CONFIGURATION AND ENVIRONMENT SETUP
# ==========================================
# Purpose: Load shared configuration and set up parallel processing environment
# Implementation: Sources configuration utilities and optimizes threading
# Input: 00_config_utils.sh configuration file
# Output: Global variables and optimized environment settings
# Edge cases: Falls back to default temp directory if TMPDIR not set

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps/construct_matrix"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# Set parallel processing parameters for optimal performance
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export TMPDIR="${TMPDIR:-/tmp}"

# ==========================================
# OPTIMIZED BEDTOOLS INTERSECT FUNCTION
# ==========================================
# Purpose: Performs high-performance intersection between peaks and fragments
# Implementation: Uses optimized sorting, temporary directories, and parallel processing
# Algorithm:
#   1. Validate input files and create optimized temporary directory
#   2. Check if fragments are pre-sorted, sort if necessary with parallel processing
#   3. Run bedtools intersect with sorted inputs for maximum performance
#   4. Validate output and clean up temporary files
# Input: peaks_file (sorted BED), fragments_file (BED), output_file (path)
# Output: Intersection results file with peak-fragment overlaps
# Edge cases: Handles compressed files, fallback sorting, memory optimization
# Returns: 0 on success, 1 on failure
run_bedtools_intersect_optimized() {
    local peaks_file="$1"
    local fragments_file="$2"
    local output_file="$3"
    
    log_info "Starting optimized bedtools intersect operation"
    log_info "Peaks file: $peaks_file"
    log_info "Fragments file: $fragments_file"
    log_info "Output file: $output_file"
    
    # ==========================================
    # INPUT VALIDATION
    # ==========================================
    # Purpose: Validates existence and accessibility of input files
    # Implementation: Checks file existence before proceeding
    # Edge cases: Returns early with error if files are missing
    
    [[ ! -f "$peaks_file" ]] && { log_error "Peaks file not found: $peaks_file"; return 1; }
    [[ ! -f "$fragments_file" ]] && { log_error "Fragments file not found: $fragments_file"; return 1; }
    
    # ==========================================
    # OPTIMIZED TEMPORARY DIRECTORY SETUP
    # ==========================================
    # Purpose: Creates temporary directory in fastest available storage
    # Implementation: Prioritizes RAM disk (/dev/shm), then TMPDIR, then fallback
    # Algorithm: Tests each location for existence and write permissions
    # Edge cases: Falls back through multiple options, ensures unique directory names
    
    local temp_intersect_dir
    for temp_path in "${TMPDIR}" "/dev/shm" "${TEMP_DIR}"; do
        if [[ -d "$temp_path" && -w "$temp_path" ]]; then
            temp_intersect_dir="${temp_path}/intersect_$$"
            break
        fi
    done
    
    mkdir -p "$temp_intersect_dir" || {
        log_error "Failed to create temporary intersect directory: $temp_intersect_dir"
        return 1
    }
    
    log_info "Using temporary directory: $temp_intersect_dir"
    
    # ==========================================
    # INTELLIGENT FRAGMENT SORTING
    # ==========================================
    # Purpose: Optimizes fragment file sorting for bedtools performance
    # Implementation: Detects pre-sorted files and applies optimized sorting when needed
    # Algorithm:
    #   1. Check filename patterns to detect pre-sorted files
    #   2. If not sorted, apply parallel sorting with optimized parameters
    #   3. Handle compressed files with appropriate decompression
    #   4. Use temporary directory for sorting operations
    # Edge cases: Handles compressed files, memory-limited sorting, cleanup on failure
    
    local fragments_sorted=false
    if [[ "${fragments_file}" == *"sorted"* ]] || [[ "${fragments_file}" == *"sorted.bed"* ]]; then
        log_info "Fragments file appears to be pre-sorted, skipping sort step"
        fragments_sorted=true
    fi
    
    local fragments_input
    if [[ "$fragments_sorted" == true ]]; then
        fragments_input="$fragments_file"
    else
        fragments_input="${temp_intersect_dir}/sorted_fragments.bed"
        log_info "Sorting fragments file with optimized parameters..."
        
        # Optimized sorting with parallel processing and memory management
        local sort_cmd
        if [[ "$fragments_file" == *.gz ]]; then
            sort_cmd="zcat '$fragments_file' | LC_ALL=C sort --parallel=$SLURM_CPUS_PER_TASK -S 24G -T '$temp_intersect_dir' -k1,1 -k2,2n"
        else
            sort_cmd="LC_ALL=C sort --parallel=$SLURM_CPUS_PER_TASK -S 24G -T '$temp_intersect_dir' -k1,1 -k2,2n '$fragments_file'"
        fi
        
        eval "$sort_cmd" > "$fragments_input" || {
            log_error "Fragment sorting failed"
            rm -rf "$temp_intersect_dir"
            return 1
        }
        
        log_info "Fragment sorting completed"
        log_info "Sorted fragments size: $(du -h "$fragments_input" | cut -f1)"
    fi
    
    # ==========================================
    # OPTIMIZED BEDTOOLS INTERSECT EXECUTION
    # ==========================================
    # Purpose: Performs intersection with maximum performance optimizations
    # Implementation: Uses sorted mode when possible, falls back gracefully
    # Algorithm:
    #   1. Attempt sorted intersect mode for optimal performance
    #   2. Fall back to standard mode if sorted mode fails
    #   3. Use -wa -wb flags to output both peak and fragment information
    # Edge cases: Handles bedtools version compatibility, ensures robust fallback
    
    log_info "Running optimized bedtools intersect..."
    
    # Use bedtools intersect with sorted input for better performance
    if bedtools intersect -sorted -a "$peaks_file" -b "$fragments_input" -wa -wb > "$output_file" 2>/dev/null; then
        log_info "Used sorted intersect mode"
    else
        log_info "Falling back to standard intersect mode"
        bedtools intersect -a "$peaks_file" -b "$fragments_input" -wa -wb > "$output_file" || {
            log_error "Bedtools intersect failed"
            rm -rf "$temp_intersect_dir"
            return 1
        }
    fi
    
    # ==========================================
    # OUTPUT VALIDATION AND CLEANUP
    # ==========================================
    # Purpose: Validates intersection results and performs cleanup
    # Implementation: Checks file existence, content, and reports statistics
    # Edge cases: Handles empty results, calculates file sizes, ensures cleanup
    
    [[ ! -f "$output_file" ]] && { log_error "Output file not created"; rm -rf "$temp_intersect_dir"; return 1; }
    
    local output_lines output_size_mb
    output_lines=$(wc -l < "$output_file")
    [[ $output_lines -eq 0 ]] && { log_error "Output file is empty"; rm -rf "$temp_intersect_dir"; return 1; }
    
    output_size_mb=$(($(stat -c%s "$output_file") / 1024 / 1024))
    log_info "Intersect completed: $output_lines lines, ${output_size_mb}MB"
    
    rm -rf "$temp_intersect_dir"
    return 0
}

# ==========================================
# MAIN EXECUTION FUNCTION
# ==========================================
# Purpose: Orchestrates the complete bedtools intersect workflow
# Implementation: Sequential execution of validation, intersection, and reporting
# Input: Uses global variables for file paths and configuration
# Output: Creates intersect output file and summary report
# Edge cases: Comprehensive error handling with early exit on failures

main() {
    log_info "=== Step 3: Optimized Bedtools Intersect Operation ==="
    log_info "Started at: $(date)"
    log_info "Using $SLURM_CPUS_PER_TASK CPU cores for parallel processing"
    
    # ==========================================
    # ENVIRONMENT INITIALIZATION
    # ==========================================
    # Purpose: Loads configuration utilities and sets up execution environment
    # Implementation: Initializes logging, conda environment, and temporary directories
    # Edge cases: Ensures all required components are properly configured
    
    init_logging
    setup_conda_env
    setup_temp_dir
    
    # ==========================================
    # DISK SPACE VALIDATION
    # ==========================================
    # Purpose: Ensures sufficient disk space for intersection output
    # Implementation: Checks for minimum 5GB available space (reduced due to optimizations)
    # Edge cases: Exits if insufficient space to prevent disk full errors
    
    check_disk_space "$OUTPUT_DIR" 5
    
    # ==========================================
    # FILE PATH DEFINITION
    # ==========================================
    # Purpose: Defines input and output file paths for intersection operation
    # Implementation: Uses global configuration variables for file locations
    # Edge cases: Ensures consistent file naming and path resolution
    
    local peaks_file="$PEAKS_FILE"
    local fragments_file="$FRAGMENTS_FILE"
    local intersect_output="${OUTPUT_DIR}/${SAMPLE}_intersect_output.txt"
    
    # ==========================================
    # INPUT FILE VALIDATION
    # ==========================================
    # Purpose: Validates existence and accessibility of required input files
    # Implementation: Checks file existence before proceeding with intersection
    # Edge cases: Exits early if any required file is missing with helpful error messages
    
    if [[ ! -f "$peaks_file" ]]; then
        log_error "Sorted peaks file not found: $peaks_file"
        log_error "Please run step 1 (validate_prerequisites.sh) first"
        exit 1
    fi
    
    if [[ ! -f "$fragments_file" ]]; then
        log_error "Fragments file not found: $fragments_file"
        exit 1
    fi
    
    # ==========================================
    # BEDTOOLS VERSION VERIFICATION
    # ==========================================
    # Purpose: Ensures bedtools is available and reports version information
    # Implementation: Checks command availability and logs version for reproducibility
    # Edge cases: Exits if bedtools is not found with installation guidance
    
    if ! command -v bedtools &> /dev/null; then
        log_error "bedtools is not available. Please install bedtools or activate the correct conda environment."
        exit 1
    fi
    
    log_info "Using bedtools version: $(bedtools --version)"
    
    # ==========================================
    # OPTIMIZED INTERSECTION EXECUTION
    # ==========================================
    # Purpose: Executes the main bedtools intersect operation with performance optimizations
    # Implementation: Calls optimized intersection function with comprehensive error handling
    # Algorithm: Uses sorted intersect mode, parallel processing, and temporary storage optimization
    # Edge cases: Exits on intersection failure with detailed error reporting and cleanup
    
    if ! run_bedtools_intersect_optimized "$peaks_file" "$fragments_file" "$intersect_output"; then
        log_error "Bedtools intersect operation failed"
        exit 1
    fi
    
    # ==========================================
    # SUMMARY STATISTICS GENERATION
    # ==========================================
    # Purpose: Generates comprehensive statistics from intersection results
    # Implementation: Uses parallel processing for efficient counting and unique value extraction
    # Algorithm: Counts total intersections, unique peaks, and unique barcodes using optimized sorting
    # Edge cases: Handles large files with memory-efficient processing and parallel operations
    
    log_info "Generating summary statistics..."
    local total_intersections unique_peaks unique_barcodes
    total_intersections=$(wc -l < "$intersect_output")
    unique_peaks=$(cut -f1-3 "$intersect_output" | LC_ALL=C sort -u --parallel=$SLURM_CPUS_PER_TASK | wc -l)
    unique_barcodes=$(cut -f7 "$intersect_output" | LC_ALL=C sort -u --parallel=$SLURM_CPUS_PER_TASK | wc -l)
    
    log_info "=== Bedtools Intersect Summary ==="
    log_info "Total intersections: $total_intersections"
    log_info "Unique peaks with intersections: $unique_peaks"
    log_info "Unique barcodes with intersections: $unique_barcodes"
    log_info "Output file: $intersect_output"
    
    # Save summary to file
    cat > "${OUTPUT_DIR}/${SAMPLE}_step3_intersect_summary.txt" << EOF
Bedtools Intersect Summary
=========================
Completed at: $(date)
Total intersections: $total_intersections
Unique peaks with intersections: $unique_peaks
Unique barcodes with intersections: $unique_barcodes
Output file: $intersect_output
EOF
    
    log_info "Step 3 completed successfully at: $(date)"
    log_info "Next step: Run 04_process_overlaps.sh"
}

# ==========================================
# SCRIPT EXECUTION ENTRY POINT
# ==========================================
# Purpose: Executes main function when script is run directly (not sourced)
# Implementation: Calls main function with all command line arguments
# Edge cases: Passes all arguments to main function for parameter flexibility

main "$@"