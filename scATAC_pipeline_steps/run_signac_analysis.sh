#!/bin/bash
#SBATCH --job-name=signac_analysis
#SBATCH --output=logs/signac_analysis.out
#SBATCH --error=logs/signac_analysis.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# ============================================================================
# SCRIPT: run_signac_analysis.sh
# PURPOSE: Run Signac scATAC-seq analysis pipeline
#
# DESCRIPTION:
# This script runs the corrected Signac analysis pipeline for scATAC-seq data.
# The pipeline includes:
# - Loading individual sample matrices (Control and Mutant)
# - Creating ChromatinAssay objects with proper genome annotations
# - Merging samples with sample-specific labels
# - Adding fragment information if available
# - Quality control and filtering
# - Dimensionality reduction and clustering
# - Downstream analysis and visualization
#
# INPUT:
# - Individual sample filtered peak-barcode matrices from previous steps
# - Fragment files (optional but recommended)
#
# OUTPUT:
# - Merged Seurat object with proper sample labels
# - Quality control plots and metrics
# - Dimensionality reduction results
# - Analysis plots and summary statistics
#
# COMPUTATIONAL REQUIREMENTS:
# - Memory: 128GB (Signac with large datasets is memory-intensive)
# - CPUs: 16 cores (parallel processing in Signac/Seurat)
# - Time: ~12 hours (depends on dataset size and analysis complexity)
#
# DEPENDENCIES:
# - R packages: Signac, Seurat, GenomeInfoDb, EnsDb.Mmusculus.v79,
#   BSgenome.Mmusculus.UCSC.mm10, JASPAR2020, TFBSTools, motifmatchr,
#   chromVAR, ggplot2, patchwork, dplyr, Matrix
# - Conda environment: signac_env (or R environment with required packages)
# ============================================================================

# Set up conda environment
# Activate specialized environment with Signac and dependencies
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate signac_env

# Enable strict error handling
# -e: exit on error, -u: exit on undefined variable, -o pipefail: exit on pipe failure
set -euo pipefail

# ============================================================================
# PIPELINE INITIALIZATION AND LOGGING
# ============================================================================
echo "========================================="
echo "Signac scATAC-seq Analysis Pipeline"
echo "Start time: $(date)"
echo "========================================="

# ============================================================================
# PATH CONFIGURATION
# ============================================================================
# Set base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data"
SCRIPT_DIR="$BASE_DIR/scATAC_pipeline_steps/Signac_3"
OUTPUT_DIR="$SCRIPT_DIR"
R_SCRIPT="$SCRIPT_DIR/signac_analysis_corrected.R"

# Create output directories
mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$SCRIPT_DIR/plots"
mkdir -p "$SCRIPT_DIR/results"

# ============================================================================
# PREREQUISITE VALIDATION
# ============================================================================
echo "Validating prerequisites..."

# Check if R script exists
if [[ ! -f "$R_SCRIPT" ]]; then
    echo "ERROR: R script not found: $R_SCRIPT"
    echo "Please ensure the Signac analysis script is available"
    exit 1
fi

# Check if input matrix directories exist
CTRL_MATRIX="$BASE_DIR/chromap_final_output/R26-Nestin-Ctrl-adult_filtered_peak_bc_matrix"
MUT_MATRIX="$BASE_DIR/chromap_final_output/R26-Nestin-Mut-adult_filtered_peak_bc_matrix"

if [[ ! -d "$CTRL_MATRIX" ]]; then
    echo "ERROR: Control matrix directory not found: $CTRL_MATRIX"
    echo "Please run the matrix creation steps first"
    exit 1
fi

if [[ ! -d "$MUT_MATRIX" ]]; then
    echo "ERROR: Mutant matrix directory not found: $MUT_MATRIX"
    echo "Please run the matrix creation steps first"
    exit 1
fi

# Check if matrix files exist
for sample_dir in "$CTRL_MATRIX" "$MUT_MATRIX"; do
    if [[ ! -f "$sample_dir/matrix.mtx.gz" ]]; then
        echo "ERROR: Matrix file not found: $sample_dir/matrix.mtx.gz"
        exit 1
    fi
    if [[ ! -f "$sample_dir/features.tsv.gz" ]]; then
        echo "ERROR: Features file not found: $sample_dir/features.tsv.gz"
        exit 1
    fi
    if [[ ! -f "$sample_dir/barcodes.tsv.gz" ]]; then
        echo "ERROR: Barcodes file not found: $sample_dir/barcodes.tsv.gz"
        exit 1
    fi
done

echo "Prerequisites validated successfully"

# ============================================================================
# SYSTEM RESOURCE MONITORING
# ============================================================================
echo "System resources:"
echo "  CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "  Memory allocated: $SLURM_MEM_PER_NODE MB"
echo "  Available memory: $(free -h | grep '^Mem:' | awk '{print $7}')"
echo "  Available disk space: $(df -h $BASE_DIR | tail -1 | awk '{print $4}')"

# ============================================================================
# R ENVIRONMENT VALIDATION
# ============================================================================
echo "Validating R environment..."

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    echo "You can install with: conda install -c conda-forge r-base"
    exit 1
fi

# Check R version
R_VERSION=$(Rscript --version 2>&1 | head -1)
echo "R version: $R_VERSION"

# Test critical R packages (quick check)
echo "Checking critical R packages..."
Rscript -e "
packages <- c('Signac', 'Seurat', 'GenomeInfoDb', 'EnsDb.Mmusculus.v79', 'BSgenome.Mmusculus.UCSC.mm10', 'Matrix')
missing <- packages[!sapply(packages, requireNamespace, quietly=TRUE)]
if(length(missing) > 0) {
    cat('ERROR: Missing packages:', paste(missing, collapse=', '), '\n')
    cat('Install with: install.packages(c(', paste(paste0(\"'\", missing, \"'\"), collapse=', '), '))\n')
    quit(status=1)
} else {
    cat('All critical packages available\n')
}
"

if [[ $? -ne 0 ]]; then
    echo "ERROR: R package validation failed"
    exit 1
fi

echo "R environment validated successfully"

# ============================================================================
# RUN SIGNAC ANALYSIS
# ============================================================================
echo "Starting Signac analysis..."
echo "Working directory: $(pwd)"
echo "Script: $R_SCRIPT"
echo "Output directory: $OUTPUT_DIR"

# Change to output directory for relative paths to work correctly
cd "$OUTPUT_DIR"

# Set environment variables for the R script
export R_MAX_NUM_DLLS=500  # Increase DLL limit for large packages
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Use all allocated CPUs

# Run the R script with detailed logging
echo "========================================="
echo "Executing R analysis script..."
echo "Start time: $(date)"
echo "========================================="

# Capture both stdout and stderr, and display progress
Rscript "$R_SCRIPT" 2>&1 | tee "logs/signac_analysis_inter.log"

# Check if analysis completed successfully
R_EXIT_CODE=${PIPESTATUS[0]}

if [[ $R_EXIT_CODE -eq 0 ]]; then
    echo "========================================="
    echo "Signac analysis completed successfully!"
    echo "End time: $(date)"
    echo "========================================="

    # Display output files
    echo "Generated files:"
    if [[ -d "plots" ]]; then
        echo "Plots directory:"
        ls -la plots/ | head -20
    fi

    if [[ -d "results" ]]; then
        echo "Results directory:"
        ls -la results/ | head -10
    fi

    echo "Log files:"
    ls -la logs/ | head -10

else
    echo "========================================="
    echo "ERROR: Signac analysis failed with exit code: $R_EXIT_CODE"
    echo "End time: $(date)"
    echo "========================================="

    echo "Recent log entries:"
    tail -50 "logs/signac_analysis_$(date +%Y%m%d)_"*.log 2>/dev/null || echo "No log files found"

    echo "Troubleshooting suggestions:"
    echo "1. Check R package installations"
    echo "2. Verify input data integrity"
    echo "3. Check memory/disk space availability"
    echo "4. Review error messages in log files"

    exit 1
fi

# ============================================================================
# POST-ANALYSIS SUMMARY
# ============================================================================
echo "Analysis Summary:"
echo "  Working directory: $OUTPUT_DIR"
echo "  Input samples: Control (R26-Nestin-Ctrl-adult), Mutant (R26-Nestin-Mut-adult)"
echo "  Log files: $OUTPUT_DIR/logs/"
echo "  Plots: $OUTPUT_DIR/plots/"
echo "  Results: $OUTPUT_DIR/results/"

echo "========================================="
echo "Signac scATAC-seq Analysis Pipeline Complete"
echo "Total runtime: $(date)"
echo "========================================="