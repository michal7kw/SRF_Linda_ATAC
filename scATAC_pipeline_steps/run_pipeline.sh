#!/bin/bash
#SBATCH --job-name=scatac_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# Master Pipeline Script: scATAC-seq Peak-Cell Matrix Generation
# This script orchestrates the execution of all pipeline steps in sequence
# with proper error handling, logging, and progress tracking
# NOTE: This script should be executed with sbatch, not bash

# ==========================================
# SLURM ENVIRONMENT CHECK
# ==========================================

# Ensure script is running under SLURM
# if [[ -z "${SLURM_JOB_ID:-}" ]]; then
#     echo "ERROR: This script must be executed with sbatch, not bash directly."
#     echo "Usage: sbatch $0"
#     exit 1
# fi

# Load shared configuration and utilities
SOURCE_DIR="$(dirname "${BASH_SOURCE[0]}")"
source "${SOURCE_DIR}/00_config_utils.sh"

# Pipeline configuration
PIPELINE_STEPS=(
    "01_validate_prerequisites.sh"
    "02_estimate_resources.sh"
    "03_bedtools_intersect.sh"
    "04_process_overlaps.sh"
    "05_create_matrix.sh"
    "06_validate_output.sh"
)

PIPELINE_NAMES=(
    "Prerequisites Validation"
    "Resource Estimation"
    "Bedtools Intersect"
    "Process Overlaps"
    "Create Matrix"
    "Output Validation"
)

# Function to display pipeline overview
show_pipeline_overview() {
    echo ""
    echo "🧬 scATAC-seq Peak-Cell Matrix Generation Pipeline 🧬"
    echo "===================================================="
    echo ""
    echo "This pipeline will execute the following steps:"
    echo ""
    
    for i in "${!PIPELINE_STEPS[@]}"; do
        local step_num=$((i + 1))
        echo "  Step $step_num: ${PIPELINE_NAMES[i]}"
        echo "           ${PIPELINE_STEPS[i]}"
        echo ""
    done
    
    echo "Configuration:"
    echo "  Sample: $SAMPLE_NAME"
    echo "  Peaks file: $PEAKS_FILE"
    echo "  Fragments file: $FRAGMENTS_FILE"
    echo "  Output directory: $OUTPUT_DIR"
    echo "  Conda environment: $CONDA_ENV"
    echo ""
}

# Function to run a single pipeline step
run_pipeline_step() {
    local step_script="$1"
    local step_name="$2"
    local step_num="$3"
    local total_steps="$4"
    
    local script_path="${SOURCE_DIR}/${step_script}"
    
    echo ""
    echo "🔄 [$step_num/$total_steps] Starting: $step_name"
    echo "================================================"
    echo "Script: $step_script"
    echo "Started at: $(date)"
    echo ""
    
    # Check if script exists
    if [[ ! -f "$script_path" ]]; then
        log_error "Pipeline step script not found: $script_path"
        return 1
    fi
    
    # Make script executable
    chmod +x "$script_path"
    
    # Run the step with sbatch
    local start_time=$(date +%s)
    
    # Submit job with sbatch and wait for completion
    local job_id=$(sbatch --parsable "$script_path")
    if [[ $? -eq 0 && -n "$job_id" ]]; then
        echo "Job submitted with ID: $job_id"
        echo "Waiting for job completion..."
        
        # Wait for job to complete
        while squeue -j "$job_id" &>/dev/null; do
            sleep 10
        done
        
        # Check job exit status
        local exit_code=$(sacct -j "$job_id" --format=ExitCode --noheader --parsable2 | head -1 | cut -d':' -f1)
        if [[ "$exit_code" == "0" ]]; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        local duration_formatted=$(format_duration $duration)
        
        echo ""
        echo "✅ [$step_num/$total_steps] Completed: $step_name"
        echo "Duration: $duration_formatted"
        echo "Completed at: $(date)"
        echo ""
        
        return 0
        else
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
            local duration_formatted=$(format_duration $duration)
            
            echo ""
            echo "❌ [$step_num/$total_steps] Failed: $step_name (Exit code: $exit_code)"
            echo "Duration: $duration_formatted"
            echo "Failed at: $(date)"
            echo "Check SLURM logs: logs/${step_script%.*}_*.out and logs/${step_script%.*}_*.err"
            echo ""
            
            return 1
        fi
    else
        echo ""
        echo "❌ [$step_num/$total_steps] Failed to submit job: $step_name"
        echo "Failed at: $(date)"
        echo ""
        
        return 1
    fi
}

# Function to format duration
format_duration() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))
    
    if [[ $hours -gt 0 ]]; then
        printf "%dh %dm %ds" $hours $minutes $secs
    elif [[ $minutes -gt 0 ]]; then
        printf "%dm %ds" $minutes $secs
    else
        printf "%ds" $secs
    fi
}

# Function to check prerequisites
check_pipeline_prerequisites() {
    log_info "Checking pipeline prerequisites"
    
    # Check if all step scripts exist
    local missing_scripts=()
    
    for script in "${PIPELINE_STEPS[@]}"; do
        local script_path="${SOURCE_DIR}/${script}"
        if [[ ! -f "$script_path" ]]; then
            missing_scripts+=("$script")
        fi
    done
    
    if [[ ${#missing_scripts[@]} -gt 0 ]]; then
        log_error "Missing pipeline step scripts: ${missing_scripts[*]}"
        return 1
    fi
    
    # Check configuration
    if [[ -z "$SAMPLE_NAME" ]]; then
        log_error "SAMPLE_NAME not configured"
        return 1
    fi
    
    if [[ -z "$PEAKS_FILE" ]]; then
        log_error "PEAKS_FILE not configured"
        return 1
    fi
    
    if [[ -z "$FRAGMENTS_FILE" ]]; then
        log_error "FRAGMENTS_FILE not configured"
        return 1
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "OUTPUT_DIR not configured"
        return 1
    fi
    
    log_info "Pipeline prerequisites check passed"
    return 0
}

# Function to create pipeline summary
create_pipeline_summary() {
    local summary_file="$1"
    local total_duration="$2"
    local failed_step="$3"
    
    local status="COMPLETED"
    if [[ -n "$failed_step" ]]; then
        status="FAILED"
    fi
    
    cat > "$summary_file" << EOF
scATAC-seq Peak-Cell Matrix Pipeline Summary
===========================================
Pipeline Status: $status
Completed at: $(date)
Total Duration: $total_duration

Configuration:
- Sample: $SAMPLE_NAME
- Peaks file: $PEAKS_FILE
- Fragments file: $FRAGMENTS_FILE
- Output directory: $OUTPUT_DIR
- Conda environment: $CONDA_ENV

Pipeline Steps:
EOF
    
    for i in "${!PIPELINE_STEPS[@]}"; do
        local step_num=$((i + 1))
        local step_status="✅ COMPLETED"
        
        if [[ -n "$failed_step" && $step_num -ge $failed_step ]]; then
            if [[ $step_num -eq $failed_step ]]; then
                step_status="❌ FAILED"
            else
                step_status="⏸️ SKIPPED"
            fi
        fi
        
        echo "Step $step_num: ${PIPELINE_NAMES[i]} - $step_status" >> "$summary_file"
    done
    
    if [[ -n "$failed_step" ]]; then
        cat >> "$summary_file" << EOF

Failure Information:
- Failed at step $failed_step: ${PIPELINE_NAMES[$((failed_step-1))]}
- Check the log files for detailed error information
- You can resume the pipeline by running the failed step manually

To resume:
  cd "${SOURCE_DIR}"
  bash "${PIPELINE_STEPS[$((failed_step-1))]}"
EOF
    else
        cat >> "$summary_file" << EOF

Success! 🎉
- All pipeline steps completed successfully
- Final matrix available in: ${OUTPUT_DIR}/filtered_peak_bc_matrix
- Final report: ${OUTPUT_DIR}/final_pipeline_report.txt

Next Steps:
1. Load the matrix in your analysis software
2. Perform quality control and filtering
3. Proceed with downstream scATAC-seq analysis
EOF
    fi
    
    log_info "Pipeline summary saved to: $summary_file"
}

# Function to handle pipeline interruption
handle_interruption() {
    echo ""
    echo "⚠️ Pipeline interrupted by user"
    echo "Current progress has been saved"
    echo "You can resume by running individual step scripts"
    echo ""
    exit 130
}

# Function to display usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help          Show this help message"
    echo "  -s, --start-step N  Start from step N (1-6)"
    echo "  -e, --end-step N    End at step N (1-6)"
    echo "  -l, --list          List all pipeline steps"
    echo "  -c, --check         Check prerequisites only"
    echo "  -y, --yes           Skip confirmation prompt"
    echo ""
    echo "Examples:"
    echo "  $0                  Run complete pipeline"
    echo "  $0 -s 3             Start from step 3"
    echo "  $0 -s 2 -e 4        Run steps 2 through 4"
    echo "  $0 -c               Check prerequisites only"
    echo ""
}

# Main execution
main() {
    # Parse command line arguments
    local start_step=1
    local end_step=${#PIPELINE_STEPS[@]}
    local check_only=false
    local skip_confirmation=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -s|--start-step)
                start_step="$2"
                shift 2
                ;;
            -e|--end-step)
                end_step="$2"
                shift 2
                ;;
            -l|--list)
                show_pipeline_overview
                exit 0
                ;;
            -c|--check)
                check_only=true
                shift
                ;;
            -y|--yes)
                skip_confirmation=true
                shift
                ;;
            *)
                echo "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
    
    # Validate step range
    if [[ $start_step -lt 1 || $start_step -gt ${#PIPELINE_STEPS[@]} ]]; then
        echo "Error: Start step must be between 1 and ${#PIPELINE_STEPS[@]}"
        exit 1
    fi
    
    if [[ $end_step -lt 1 || $end_step -gt ${#PIPELINE_STEPS[@]} ]]; then
        echo "Error: End step must be between 1 and ${#PIPELINE_STEPS[@]}"
        exit 1
    fi
    
    if [[ $start_step -gt $end_step ]]; then
        echo "Error: Start step cannot be greater than end step"
        exit 1
    fi
    
    # Initialize environment
    setup_logging
    
    log_info "=== scATAC-seq Peak-Cell Matrix Pipeline ==="
    log_info "Started at: $(date)"
    
    # Set up signal handlers
    trap handle_interruption SIGINT SIGTERM
    
    # Check prerequisites
    if ! check_pipeline_prerequisites; then
        log_error "Prerequisites check failed"
        exit 1
    fi
    
    if [[ $check_only == true ]]; then
        echo "✅ Prerequisites check passed"
        exit 0
    fi
    
    # Show pipeline overview
    show_pipeline_overview
    
    # Confirmation prompt
    if [[ $skip_confirmation == false ]]; then
        echo "Do you want to proceed with the pipeline? (y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "Pipeline cancelled by user"
            exit 0
        fi
        echo ""
    fi
    
    # Initialize environment for pipeline
    activate_conda_env
    setup_temp_dir
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR" || {
        log_error "Failed to create output directory: $OUTPUT_DIR"
        exit 1
    }
    
    # Run pipeline steps
    local pipeline_start_time=$(date +%s)
    local failed_step=""
    
    for i in $(seq $((start_step-1)) $((end_step-1))); do
        local step_num=$((i + 1))
        local step_script="${PIPELINE_STEPS[i]}"
        local step_name="${PIPELINE_NAMES[i]}"
        
        if ! run_pipeline_step "$step_script" "$step_name" "$step_num" "$end_step"; then
            failed_step="$step_num"
            break
        fi
    done
    
    local pipeline_end_time=$(date +%s)
    local total_duration=$((pipeline_end_time - pipeline_start_time))
    local total_duration_formatted=$(format_duration $total_duration)
    
    # Create pipeline summary
    local summary_file="${OUTPUT_DIR}/pipeline_summary.txt"
    create_pipeline_summary "$summary_file" "$total_duration_formatted" "$failed_step"
    
    # Final status
    echo ""
    echo "================================================"
    
    if [[ -n "$failed_step" ]]; then
        echo "❌ Pipeline FAILED at step $failed_step"
        echo "Total duration: $total_duration_formatted"
        echo "Summary: $summary_file"
        echo ""
        echo "To resume, run the failed step manually:"
        echo "  bash \"${SOURCE_DIR}/${PIPELINE_STEPS[$((failed_step-1))]}\""
        echo ""
        exit 1
    else
        echo "🎉 Pipeline COMPLETED successfully!"
        echo "Total duration: $total_duration_formatted"
        echo "Summary: $summary_file"
        echo "Final matrix: ${OUTPUT_DIR}/filtered_peak_bc_matrix"
        echo ""
        echo "Your scATAC-seq peak-cell matrix is ready for analysis! 🧬"
        echo ""
        exit 0
    fi
}

# Run main function with all arguments
main "$@"