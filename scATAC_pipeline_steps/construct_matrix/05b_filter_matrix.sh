#!/bin/bash
#SBATCH --job-name=05b_filter_matrix
#SBATCH --output=logs/05b_filter_matrix_%a.out
#SBATCH --error=logs/05b_filter_matrix_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=workq

# Step 5b: Filter Peak-Barcode Matrix
# This script applies quality control filtering to the raw peak-barcode matrix
# and creates a properly filtered dataset for downstream analysis

# ==========================================
# SLURM ENVIRONMENT CHECK
# ==========================================

# Ensure script is running under SLURM
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: This script must be executed with sbatch, not bash directly."
    echo "Usage: sbatch $0"
    exit 1
fi

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps/construct_matrix"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# ==========================================
# SAMPLE CONFIGURATION
# ==========================================

SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "========================================="
echo "Step 5b: Matrix Filtering for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# ==========================================
# PATH SETUP
# ==========================================

# Input: Raw matrix from step 5
INPUT_MATRIX_DIR="${OUTPUT_DIR}/${SAMPLE}_filtered_peak_bc_matrix"

# Output: Filtered matrix in separate directory
FILTERED_MATRIX_DIR="${OUTPUT_DIR}/${SAMPLE}_filtered_peak_bc_matrix_qc"

# Create output directory
mkdir -p "$FILTERED_MATRIX_DIR"

# ==========================================
# PREREQUISITE VALIDATION
# ==========================================

# Check if raw matrix exists
if [[ ! -d "$INPUT_MATRIX_DIR" ]]; then
    log_error "Raw matrix directory not found: $INPUT_MATRIX_DIR"
    log_error "Please run step 5 (05_create_matrix.sh) first"
    exit 1
fi

if [[ ! -f "$INPUT_MATRIX_DIR/matrix.mtx.gz" ]]; then
    log_error "Matrix file not found: $INPUT_MATRIX_DIR/matrix.mtx.gz"
    exit 1
fi

log_info "Prerequisites verified - raw matrix files found"

# ==========================================
# R SCRIPT FOR MATRIX FILTERING
# ==========================================

cat > "${FILTERED_MATRIX_DIR}/filter_matrix.R" << 'EOF'
#!/usr/bin/env Rscript

# ============================================================================
# MATRIX FILTERING SCRIPT FOR scATAC-seq DATA
# ============================================================================
# PURPOSE: Apply quality control filtering to peak-barcode matrix
# INPUT: Raw 10X format matrix (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
# OUTPUT: Filtered matrix with same format but high-quality cells/peaks only
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(Matrix)        # Sparse matrix operations
})

# Get sample name from environment
sample_name <- Sys.getenv("SAMPLE", "unknown")
cat("Filtering matrix for sample:", sample_name, "\n")

# ============================================================================
# PATH CONFIGURATION
# ============================================================================

input_dir <- paste0("../", sample_name, "_filtered_peak_bc_matrix")
output_dir <- "."

# ============================================================================
# LOAD RAW MATRIX DATA
# ============================================================================

cat("Loading raw peak-barcode matrix...\n")

# Load 10X format files
matrix_path <- file.path(input_dir, "matrix.mtx.gz")
features_path <- file.path(input_dir, "features.tsv.gz")
barcodes_path <- file.path(input_dir, "barcodes.tsv.gz")

# Validate input files
if(!file.exists(matrix_path)) {
    stop("Matrix file not found: ", matrix_path)
}
if(!file.exists(features_path)) {
    stop("Features file not found: ", features_path)
}
if(!file.exists(barcodes_path)) {
    stop("Barcodes file not found: ", barcodes_path)
}

# Read matrix and metadata
peak_matrix <- readMM(matrix_path)
features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(peak_matrix) <- features[,1]  # Peak coordinates
colnames(peak_matrix) <- barcodes[,1]  # Cell barcodes

cat("Raw matrix dimensions:", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# ============================================================================
# CALCULATE QC METRICS
# ============================================================================

cat("Calculating QC metrics...\n")

# Calculate peaks per cell and cells per peak
cell_peak_counts <- Matrix::colSums(peak_matrix > 0)
peak_cell_counts <- Matrix::rowSums(peak_matrix > 0)
cell_total_counts <- Matrix::colSums(peak_matrix)

# ============================================================================
# QUALITY CONTROL STATISTICS
# ============================================================================

cat("=== RAW DATA STATISTICS ===\n")
cat("Peak count distribution per cell:\n")
cat("  Min:", min(cell_peak_counts), "\n")
cat("  Max:", max(cell_peak_counts), "\n") 
cat("  Median:", median(cell_peak_counts), "\n")
cat("  Mean:", round(mean(cell_peak_counts), 2), "\n")
cat("  Quantiles (10%, 25%, 50%, 75%, 90%, 95%, 98%):\n")
print(quantile(cell_peak_counts, c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98)))

cat("Cell distribution by peak count:\n")
cat("  Cells with >0 peaks:", sum(cell_peak_counts > 0), "\n")
cat("  Cells with >10 peaks:", sum(cell_peak_counts > 10), "\n")
cat("  Cells with >100 peaks:", sum(cell_peak_counts > 100), "\n")
cat("  Cells with >500 peaks:", sum(cell_peak_counts > 500), "\n")
cat("  Cells with >1000 peaks:", sum(cell_peak_counts > 1000), "\n")

cat("Peak distribution:\n")
cat("  Peaks in >0 cells:", sum(peak_cell_counts > 0), "\n")
cat("  Peaks in >10 cells:", sum(peak_cell_counts > 10), "\n")
cat("  Peaks in >100 cells:", sum(peak_cell_counts > 100), "\n")
cat("  Peaks in >1000 cells:", sum(peak_cell_counts > 1000), "\n")

# ============================================================================
# DEFINE FILTERING THRESHOLDS
# ============================================================================

cat("=== CALCULATING FILTERING THRESHOLDS ===\n")

# Cell filtering: Focus on high-quality cells
# Minimum peaks per cell: Use 100 as baseline, adapt if needed
min_peaks_per_cell <- max(100, quantile(cell_peak_counts[cell_peak_counts > 0], 0.05))

# Maximum peaks per cell: Remove potential doublets (98th percentile)
max_peaks_per_cell <- quantile(cell_peak_counts, 0.98)

# Peak filtering: Remove very rare peaks that likely represent noise
# Minimum cells per peak: At least 0.1% of cells or 10 cells minimum
min_cells_per_peak <- max(10, ceiling(ncol(peak_matrix) * 0.001))

cat("Filtering thresholds:\n")
cat("  Min peaks per cell:", min_peaks_per_cell, "\n")
cat("  Max peaks per cell:", max_peaks_per_cell, "\n")
cat("  Min cells per peak:", min_cells_per_peak, "\n")

# ============================================================================
# APPLY FILTERS
# ============================================================================

cat("=== APPLYING FILTERS ===\n")

# Identify cells and peaks that pass QC
good_cells <- which(cell_peak_counts >= min_peaks_per_cell & 
                   cell_peak_counts <= max_peaks_per_cell)
good_peaks <- which(peak_cell_counts >= min_cells_per_peak)

cat("Cells passing filter:", length(good_cells), "/", ncol(peak_matrix), 
    "(", round(100*length(good_cells)/ncol(peak_matrix), 1), "%)\n")
cat("Peaks passing filter:", length(good_peaks), "/", nrow(peak_matrix), 
    "(", round(100*length(good_peaks)/nrow(peak_matrix), 1), "%)\n")

# Create filtered matrix
filtered_matrix <- peak_matrix[good_peaks, good_cells]
filtered_features <- features[good_peaks, ]
filtered_barcodes <- barcodes[good_cells, ]

cat("Filtered matrix dimensions:", nrow(filtered_matrix), "peaks x", ncol(filtered_matrix), "cells\n")

# ============================================================================
# FILTERED DATA STATISTICS
# ============================================================================

cat("=== FILTERED DATA STATISTICS ===\n")
filtered_cell_peaks <- Matrix::colSums(filtered_matrix > 0)
filtered_cell_counts <- Matrix::colSums(filtered_matrix)

cat("Filtered peak count distribution per cell:\n")
cat("  Min:", min(filtered_cell_peaks), "\n")
cat("  Max:", max(filtered_cell_peaks), "\n")
cat("  Median:", median(filtered_cell_peaks), "\n")
cat("  Mean:", round(mean(filtered_cell_peaks), 2), "\n")
cat("  Quantiles (25%, 50%, 75%, 95%):\n")
print(quantile(filtered_cell_peaks, c(0.25, 0.50, 0.75, 0.95)))

# ============================================================================
# SAVE FILTERED MATRIX IN 10X FORMAT
# ============================================================================

cat("=== SAVING FILTERED MATRIX ===\n")

# Save matrix in Matrix Market format
Matrix::writeMM(filtered_matrix, file.path(output_dir, "matrix.mtx"))

# Save features (peaks)
write.table(filtered_features, 
           file.path(output_dir, "features.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save barcodes (cells)  
write.table(filtered_barcodes,
           file.path(output_dir, "barcodes.tsv"), 
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Compress files to match 10X format
system("gzip -f matrix.mtx")
system("gzip -f features.tsv") 
system("gzip -f barcodes.tsv")

# ============================================================================
# CREATE FILTERING SUMMARY
# ============================================================================

# Create summary statistics
summary_stats <- data.frame(
    Metric = c("Raw_peaks", "Raw_cells", "Filtered_peaks", "Filtered_cells",
               "Cells_removed", "Peaks_removed", "Cells_retained_pct", "Peaks_retained_pct",
               "Min_peaks_per_cell", "Max_peaks_per_cell", "Min_cells_per_peak"),
    Value = c(nrow(peak_matrix), ncol(peak_matrix), 
              nrow(filtered_matrix), ncol(filtered_matrix),
              ncol(peak_matrix) - ncol(filtered_matrix),
              nrow(peak_matrix) - nrow(filtered_matrix),
              round(100*ncol(filtered_matrix)/ncol(peak_matrix), 1),
              round(100*nrow(filtered_matrix)/nrow(peak_matrix), 1),
              min_peaks_per_cell, max_peaks_per_cell, min_cells_per_peak)
)

write.table(summary_stats,
           file.path(output_dir, paste0(sample_name, "_filtering_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Create detailed QC metrics for each cell
cell_qc <- data.frame(
    Barcode = colnames(filtered_matrix),
    n_peaks = Matrix::colSums(filtered_matrix > 0),
    total_counts = Matrix::colSums(filtered_matrix),
    log10_n_peaks = log10(Matrix::colSums(filtered_matrix > 0) + 1),
    log10_total_counts = log10(Matrix::colSums(filtered_matrix) + 1)
)

write.table(cell_qc,
           file.path(output_dir, paste0(sample_name, "_cell_qc_metrics.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("=== FILTERING COMPLETE ===\n")
cat("Files created:\n")
cat("  - matrix.mtx.gz: Filtered sparse matrix\n")
cat("  - features.tsv.gz: Filtered peak coordinates\n")
cat("  - barcodes.tsv.gz: Filtered cell barcodes\n")
cat("  -", paste0(sample_name, "_filtering_summary.tsv: Filtering statistics\n"))
cat("  -", paste0(sample_name, "_cell_qc_metrics.tsv: Per-cell QC metrics\n"))

cat("Matrix size reduction:\n")
cat("  Cells:", ncol(peak_matrix), "->", ncol(filtered_matrix), 
    "(", round(100*ncol(filtered_matrix)/ncol(peak_matrix), 1), "% retained)\n")
cat("  Peaks:", nrow(peak_matrix), "->", nrow(filtered_matrix),
    "(", round(100*nrow(filtered_matrix)/nrow(peak_matrix), 1), "% retained)\n")

EOF

# ==========================================
# EXECUTE FILTERING
# ==========================================

log_info "Running matrix filtering analysis..."
cd "$FILTERED_MATRIX_DIR"
export SAMPLE="$SAMPLE"

# Activate conda environment with R
log_info "Activating peak_calling_new conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate peak_calling_new

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    log_error "Rscript not found. Please install R in peak_calling_new environment"
    exit 1
fi

# Run the filtering script
log_info "Starting R filtering analysis for $SAMPLE..."
Rscript "filter_matrix.R"

if [[ $? -eq 0 ]]; then
    log_info "Matrix filtering completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="${SAMPLE}_filtering_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo ""
        echo "=== FILTERING SUMMARY ==="
        cat "$SUMMARY_FILE"
        echo ""
    fi
    
    # Validate output files
    if [[ -f "matrix.mtx.gz" && -f "features.tsv.gz" && -f "barcodes.tsv.gz" ]]; then
        log_info "All output files created successfully"
        
        # Show file sizes
        echo "Output file sizes:"
        ls -lh matrix.mtx.gz features.tsv.gz barcodes.tsv.gz
        
    else
        log_error "Some output files are missing"
        exit 1
    fi
    
else
    log_error "Matrix filtering failed. Check the R script and dependencies."
    exit 1
fi

echo ""
echo "========================================="
echo "Step 5b complete for $SAMPLE"
echo "Matrix filtering finished"
echo "End time: $(date)"
echo "========================================="
echo ""
echo "Filtered matrix saved to: $FILTERED_MATRIX_DIR"
echo "Files created:"
echo "  - matrix.mtx.gz: Filtered sparse matrix"
echo "  - features.tsv.gz: Filtered peak coordinates" 
echo "  - barcodes.tsv.gz: Filtered cell barcodes"
echo "  - ${SAMPLE}_filtering_summary.tsv: Filtering statistics"
echo "  - ${SAMPLE}_cell_qc_metrics.tsv: Per-cell QC metrics"
echo ""
echo "To use filtered data in downstream analysis, update input paths to:"
echo "  $FILTERED_MATRIX_DIR"