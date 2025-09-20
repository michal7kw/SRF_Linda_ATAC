#!/bin/bash

# =============================================================================
# scATAC-seq Pipeline Configuration and Utilities
# =============================================================================
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# Description: Shared configuration, logging, and utility functions for peak-cell matrix pipeline
#
# PURPOSE:
# This script provides centralized configuration management and utility functions
# for the scATAC-seq peak-cell matrix construction pipeline. It establishes
# global variables, logging infrastructure, and common functions used across
# all pipeline scripts.
#
# KEY FEATURES:
# - Sample configuration and path management
# - Comprehensive logging system with multiple log levels
# - Resource management and validation functions
# - Temporary directory management with automatic cleanup
# - Tool validation and environment setup
# - Error handling and progress tracking
# =============================================================================

# BASH CONFIGURATION:
# Enable strict error handling to catch issues early
# -e: Exit on any command failure
# -u: Exit on undefined variable usage
# -o pipefail: Exit on any failure in a pipeline
set -euo pipefail

# =============================================================================
# GLOBAL CONFIGURATION SECTION
# =============================================================================
# PURPOSE: Define pipeline-wide configuration parameters and sample selection
# INPUT: SLURM_ARRAY_TASK_ID environment variable or command line argument
# OUTPUT: Sets SAMPLE variable for use throughout pipeline
# LOGIC: Supports both SLURM array jobs and manual execution modes
# =============================================================================

# SAMPLE CONFIGURATION:
# Define all samples to be processed in the pipeline
# Each sample corresponds to a different experimental condition
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")

# SAMPLE SELECTION LOGIC:
# Determine which sample to process based on execution context
# Priority: 1) SLURM array task ID, 2) Command line argument, 3) Error
# This allows flexible execution in both batch and interactive modes
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array job: Use task ID as index into SAMPLES array
    SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
elif [[ $# -gt 0 ]]; then
    # Manual execution: Use first command line argument as sample name
    SAMPLE="$1"
else
    # ERROR CASE: No sample specified - provide helpful error message
    echo "ERROR: No sample specified. Use SLURM_ARRAY_TASK_ID or provide sample name as argument."
    echo "Available samples: ${SAMPLES[*]}"
    exit 1
fi

# DIRECTORY PATH CONFIGURATION:
# Establish all major directory paths used throughout the pipeline
# These paths are designed for the BeeGFS scratch filesystem structure
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data"
OUTPUT_DIR="$BASE_DIR/chromap_final_output"  # Main output directory for results
SCRIPT_DIR="$BASE_DIR/scATAC_pipeline_steps"  # Directory containing pipeline scripts

# TEMPORARY DIRECTORY MANAGEMENT:
# Create unique temporary directories to avoid conflicts between concurrent jobs
# Uses BeeGFS scratch space instead of /tmp for better performance and capacity
# LOGIC: Generate unique temp dir based on SLURM job context or timestamp
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # SLURM array job: Include both job ID and array task ID for uniqueness
    TEMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/temp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # SLURM single job: Use job ID only
    TEMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/temp/${SLURM_JOB_ID}"
else
    # Manual execution: Use timestamp to ensure uniqueness
    TEMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/temp/manual_$(date +%s)"
fi

# SAMPLE-SPECIFIC FILE PATH CONFIGURATION:
# Define all input and output file paths based on the selected sample
# These paths follow the standard scATAC-seq pipeline output structure
GENOME_SIZES="$OUTPUT_DIR/qc/genome.chrom.sizes"                    # Reference genome chromosome sizes
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"           # Sorted peak regions for the sample
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"                     # Processed reads in BED format
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"   # Fragment file from chromap
MATRIX_DIR="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"                   # Output directory for peak-barcode matrix

# RESOURCE CONFIGURATION:
# Define default resource limits and requirements for pipeline operations
# These values are optimized for typical scATAC-seq datasets
DEFAULT_SORT_MEM="32G"          # Default memory allocation for sort operations
MIN_DISK_SPACE_GB=300           # Minimum required disk space (GB) for pipeline execution
MIN_PROCESSING_SPACE_GB=200     # Minimum processing space (GB) for intermediate files

# =============================================================================
# LOGGING SYSTEM
# =============================================================================
# PURPOSE: Provide comprehensive logging infrastructure for pipeline monitoring
# FEATURES: Multiple log levels, timestamping, file and console output
# OUTPUT: Structured log messages with timestamps and severity levels
# USAGE: Called by all pipeline scripts for consistent logging
# =============================================================================

# LOGGING INITIALIZATION:
# Set up logging infrastructure and create log file for the current sample
# INPUT: Uses global SAMPLE and SCRIPT_DIR variables
# OUTPUT: Creates log directory and initializes LOG_FILE variable
# SIDE EFFECTS: Creates/clears log file, sets up logging for current session
init_logging() {
    local log_dir="$SCRIPT_DIR/logs"
    mkdir -p "$log_dir"
    
    # Set log file path with sample-specific naming
    LOG_FILE="$log_dir/06b_${SAMPLE}_progress.log"
    
    # Create or clear log file for fresh start
    > "$LOG_FILE"
    
    logger "Logging initialized for sample: $SAMPLE"
    logger "Log file: $LOG_FILE"
}

# CORE LOGGING FUNCTIONS:
# Provide different log levels for various types of messages
# All functions use tee to write to both console and log file simultaneously
# Fallback to stderr if LOG_FILE is not set

# INFORMATIONAL LOGGING:
# Standard information messages for normal pipeline operations
log_info() {
    local message="$1"
    local timestamp="$(date '+%H:%M:%S')"
    echo "[INFO $timestamp] $message" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# GENERAL LOGGER:
# Alias for log_info - maintains backward compatibility
logger() {
    local message="$1"
    local timestamp="$(date '+%H:%M:%S')"
    echo "[INFO $timestamp] $message" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# ERROR LOGGING:
# Critical error messages that indicate pipeline failures
log_error() {
    local message="$1"
    local timestamp="$(date '+%H:%M:%S')"
    echo "[ERROR $timestamp] $message" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# WARNING LOGGING:
# Non-critical issues that may affect pipeline performance or results
log_warning() {
    local message="$1"
    local timestamp="$(date '+%H:%M:%S')"
    echo "[WARNING $timestamp] $message" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# PROGRESS TRACKING:
# Track completion percentage for long-running operations
# INPUT: stage name and completion percentage
log_progress() {
    local stage="$1"
    local percent="$2"
    local timestamp="$(date '+%H:%M:%S')"
    echo "[PROGRESS $timestamp] [$stage] $percent% complete" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# STEP TRACKING:
# Track major pipeline steps with their status (START, COMPLETE, FAILED)
# Useful for monitoring pipeline execution flow
log_step() {
    local step_name="$1"
    local status="$2"  # START, COMPLETE, FAILED
    local timestamp="$(date '+%H:%M:%S')"
    echo "[STEP $timestamp] [$step_name] $status" | tee -a "${LOG_FILE:-/dev/stderr}"
}

# RESOURCE USAGE REPORTING:
# Monitor and log current system resource utilization
# PURPOSE: Help diagnose performance issues and resource constraints
# OUTPUT: Memory usage and temporary directory disk usage
# EDGE CASES: Handles missing tools and directories gracefully
report_resource_usage() {
    # Check memory usage if 'free' command is available
    if command -v free &> /dev/null; then
        logger "Memory usage: $(free -h | grep '^Mem:' | awk '{print $3"/"$2}')"
    fi
    # Check temporary directory disk usage if it exists
    if [[ -d "$TEMP_DIR" ]]; then
        logger "Disk usage (temp): $(du -sh $TEMP_DIR 2>/dev/null | cut -f1 || echo 'N/A')"
    fi
}

# =============================================================================
# UTILITY FUNCTIONS SECTION
# =============================================================================
# PURPOSE: Provide common utility functions for file management and system setup
# FEATURES: Temporary directory management, disk space checking, environment setup
# USAGE: Called by pipeline scripts for consistent resource management
# =============================================================================

# TEMPORARY DIRECTORY SETUP:
# Create and configure temporary directory for pipeline operations
# PURPOSE: Provide isolated workspace for intermediate files with automatic cleanup
# INPUT: Uses global TEMP_DIR variable
# OUTPUT: Creates directory structure and sets environment variables
# SIDE EFFECTS: Sets up EXIT trap for cleanup, exports TMPDIR
# EDGE CASES: Handles existing directories and environment variables
setup_temp_dir() {
    logger "Setting up temporary directory: $TEMP_DIR"
    
    # Create temporary directory with appropriate permissions
    mkdir -p "$TEMP_DIR"
    chmod 755 "$TEMP_DIR"
    
    # Clean up any existing temporary files from previous runs
    # This prevents conflicts and ensures clean state
    rm -rf "$TEMP_DIR"/*
    
    # Set up automatic cleanup on script exit (success or failure)
    trap "cleanup_temp_dir" EXIT
    
    # Configure system TMPDIR environment variable for tools that use it
    if [[ -z "${TMPDIR:-}" ]]; then
        export TMPDIR="$TEMP_DIR"
        logger "Set TMPDIR to: $TMPDIR"
    else
        logger "Using existing TMPDIR: $TMPDIR"
    fi
    
    logger "Temporary directory setup complete"
}

# TEMPORARY DIRECTORY CLEANUP:
# Remove temporary directory and all contents
# PURPOSE: Clean up resources on script exit to prevent disk space issues
# INPUT: Uses global TEMP_DIR variable
# EDGE CASES: Safely handles non-existent directories
cleanup_temp_dir() {
    if [[ -d "$TEMP_DIR" ]]; then
        logger "Cleaning up temporary directory: $TEMP_DIR"
        rm -rf "$TEMP_DIR"
    fi
}

# DISK SPACE VALIDATION:
# Check if sufficient disk space is available for pipeline operations
# PURPOSE: Prevent pipeline failures due to insufficient disk space
# INPUT: directory path and required space in GB
# OUTPUT: Returns 0 if sufficient space, 1 if insufficient
# ALGORITHM: Uses df command to get available space, converts to GB
# EDGE CASES: Handles directories that don't exist or are inaccessible
check_disk_space() {
    local dir="$1"
    local required_gb="$2"
    local available_kb=$(df -k "$dir" | awk 'NR==2 {print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [[ $available_gb -lt $required_gb ]]; then
        log_error "Insufficient disk space in $dir: ${available_gb}GB available, ${required_gb}GB required"
        return 1
    fi
    
    logger "Disk space check passed: ${available_gb}GB available in $dir"
    return 0
}

# CONDA ENVIRONMENT SETUP:
# Activate the required conda environment for pipeline tools
# PURPOSE: Ensure all required bioinformatics tools are available
# INPUT: None (uses hardcoded paths specific to the cluster environment)
# OUTPUT: Activates conda environment with required tools
# DEPENDENCIES: Requires conda installation at specified path
# EDGE CASES: Handles missing conda installation gracefully
setup_conda_env() {
    logger "Setting up conda environment..."
    
    # Check if conda activation script exists at expected location
    if [[ -f "/beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin//activate" ]]; then
        # Source conda and activate the environment with bioinformatics tools
        source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/bin/activate
        conda activate peak_calling_new
        logger "Conda environment 'peak_calling_new' activated"
    else
        log_error "Conda activation script not found"
        return 1
    fi
}

# =============================================================================
# VALIDATION FUNCTIONS SECTION
# =============================================================================
# PURPOSE: Provide comprehensive validation for files, tools, and system state
# FEATURES: File existence/size validation, tool availability checking
# OUTPUT: Returns 0 for success, 1 for failure with detailed error messages
# USAGE: Called before pipeline operations to ensure prerequisites are met
# =============================================================================

# FILE VALIDATION:
# Verify file exists and meets minimum size requirements
# PURPOSE: Ensure input files are present and not corrupted/empty
# INPUT: file path, minimum size (optional, default 1), description (optional)
# OUTPUT: Returns 0 if valid, 1 if invalid
# ALGORITHM: Checks existence with -f test, gets size with stat command
# EDGE CASES: Handles missing files, inaccessible files, stat command failures
validate_file() {
    local file="$1"
    local min_size="${2:-1}"        # Default minimum size: 1 byte
    local description="${3:-"File"}" # Default description: "File"
    
    # Check if file exists and is a regular file
    if [[ ! -f "$file" ]]; then
        log_error "$description not found: $file"
        return 1
    fi
    
    # Check file size meets minimum requirement
    local size=$(stat -c%s "$file" 2>/dev/null || echo 0)
    if [[ $size -lt $min_size ]]; then
        log_error "$description is too small ($size bytes): $file"
        return 1
    fi
    
    logger "$description validated: $file ($size bytes)"
    return 0
}

# TOOL AVAILABILITY VALIDATION:
# Check that all required command-line tools are available in PATH
# PURPOSE: Prevent pipeline failures due to missing dependencies
# INPUT: None (uses predefined list of required tools)
# OUTPUT: Returns 0 if all tools available, 1 if any missing
# ALGORITHM: Uses 'command -v' to check each tool's availability
# EDGE CASES: Handles tools not in PATH, permission issues
validate_tools() {
    logger "Validating required tools..."
    
    # Define list of essential tools required by the pipeline
    local tools=("bedtools" "sort" "awk" "gzip" "zcat")
    local missing_tools=()
    
    # Check each tool's availability using command -v
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    # Report results
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        return 1
    fi
    
    logger "All required tools are available"
    return 0
}

# DYNAMIC MEMORY ALLOCATION:
# Calculate optimal memory allocation for sort operations based on system resources
# PURPOSE: Optimize performance while preventing out-of-memory errors
# INPUT: None (queries system memory automatically)
# OUTPUT: Memory specification string (e.g., "32G", "16G")
# ALGORITHM: Uses 70% of available memory, caps at default maximum
# EDGE CASES: Handles systems without 'free' command, low-memory systems
get_sort_memory() {
    local available_mem
    
    # Try to get available memory using 'free' command
    if command -v free &> /dev/null; then
        # Calculate 70% of available memory to leave room for other processes
        available_mem=$(free -g | awk '/^Mem:/ {print int($7 * 0.7)}')
    else
        # Fallback for systems without 'free' command
        available_mem=16
    fi
    
    local sort_mem="$DEFAULT_SORT_MEM"
    
    # Adjust memory allocation for low-memory systems
    if [[ $available_mem -lt 32 ]]; then
        sort_mem="${available_mem}G"
        logger "Adjusting sort memory to ${sort_mem} based on available system memory"
    else
        logger "Using default sort memory: ${sort_mem}"
    fi
    
    echo "$sort_mem"
}

# =============================================================================
# INITIALIZATION AND EXPORTS
# =============================================================================
# PURPOSE: Make all configuration and functions available to other pipeline scripts
# MECHANISM: Export variables and functions to environment
# USAGE: Other scripts source this file to inherit all configuration
# =============================================================================

# EXPORT CONFIGURATION VARIABLES:
# Make all path and sample configuration available to child scripts
export SAMPLE BASE_DIR OUTPUT_DIR SCRIPT_DIR TEMP_DIR
export GENOME_SIZES PEAKS_FILE READS_FILE FRAGMENTS_FILE MATRIX_DIR
export LOG_FILE

# EXPORT LOGGING FUNCTIONS:
# Make all logging functions available for consistent logging across scripts
export -f logger log_error log_warning log_progress log_step report_resource_usage

# EXPORT UTILITY FUNCTIONS:
# Make all utility and validation functions available to other scripts
export -f setup_temp_dir cleanup_temp_dir check_disk_space setup_conda_env
export -f validate_file validate_tools get_sort_memory

# INITIALIZATION COMPLETE:
# Log successful loading of configuration for the current sample
logger "Configuration and utilities loaded for sample: $SAMPLE"