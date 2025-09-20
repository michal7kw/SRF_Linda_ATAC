#!/bin/bash
#SBATCH --job-name=dim_reduction
#SBATCH --output=logs/09_dim_reduction_%a.out
#SBATCH --error=logs/09_dim_reduction_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=workq

# ============================================================================
# SCRIPT: 09_dimensionality_reduction.sh
# PURPOSE: Perform dimensionality reduction and clustering on scATAC-seq data
# 
# DESCRIPTION:
# This script implements a comprehensive dimensionality reduction pipeline for
# single-cell ATAC-seq data, including:
# - TF-IDF normalization for sparse chromatin accessibility data
# - Latent Semantic Indexing (LSI) for initial dimensionality reduction
# - Principal Component Analysis (PCA) on LSI components
# - UMAP embedding for visualization
# - K-means clustering for cell type identification
# - Quality assessment using silhouette analysis
# 
# INPUT:
# - Filtered peak-barcode matrix from step 5 (matrix.mtx.gz format)
# - Features and barcodes files (features.tsv.gz, barcodes.tsv.gz)
# 
# OUTPUT:
# - LSI embeddings (reduced dimensional representation)
# - PCA embeddings (further dimensionality reduction)
# - UMAP coordinates (2D visualization)
# - Cluster assignments with silhouette scores
# - Visualization plots (UMAP clusters, PCA variance)
# - Summary statistics and filtered matrix
# 
# COMPUTATIONAL REQUIREMENTS:
# - Memory: 64GB (for large sparse matrix operations)
# - CPUs: 16 cores (parallel processing in R libraries)
# - Time: ~8 hours (depends on matrix size and complexity)
# 
# DEPENDENCIES:
# - R packages: Matrix, irlba, uwot, cluster, ggplot2, dplyr, RColorBrewer
# - Conda environment: dm_reduce
# ============================================================================

# Set up conda environment
# Activate specialized environment with dimensionality reduction tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate dm_reduce

# Enable strict error handling
# -e: exit on error, -u: exit on undefined variable, -o pipefail: exit on pipe failure
set -euo pipefail

# ============================================================================
# SAMPLE CONFIGURATION AND PATH SETUP
# ============================================================================
# Define sample array for batch processing
# Each array job processes one sample independently
# Array indices: 0=Control, 1=Mutant
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# Set base output directory
# Contains all processed scATAC-seq data from previous pipeline steps
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

# ============================================================================
# PIPELINE INITIALIZATION AND LOGGING
# ============================================================================
echo "========================================="
echo "Step 9: Dimensionality Reduction for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# ============================================================================
# OUTPUT DIRECTORY PREPARATION
# ============================================================================
# Create dedicated directory for dimensionality reduction results
# This will contain LSI/PCA embeddings, UMAP coordinates, clusters, and plots
mkdir -p "$OUTPUT_DIR/dimensionality_reduction"

# ============================================================================
# PREREQUISITE VALIDATION
# ============================================================================
# Verify that required input files exist from previous pipeline steps
# The peak-barcode matrix is essential for dimensionality reduction analysis

# Check for peak-barcode matrix directory from step 5
# This directory should contain matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
PEAK_BC_MATRIX="$OUTPUT_DIR/${SAMPLE}_filtered_peak_bc_matrix_qc"
if [[ ! -d "$PEAK_BC_MATRIX" ]]; then
    echo "ERROR: Peak-barcode matrix directory not found: $PEAK_BC_MATRIX"
    echo "Please run step 5 (05_create_matrix.sh) first"
    echo "This step requires the filtered peak-barcode matrix for analysis"
    exit 1
fi

# Verify the main matrix file exists
# matrix.mtx.gz contains the sparse peak accessibility counts per cell
if [[ ! -f "$PEAK_BC_MATRIX/matrix.mtx.gz" ]]; then
    echo "ERROR: Matrix file not found: $PEAK_BC_MATRIX/matrix.mtx.gz"
    echo "The sparse matrix file is required for dimensionality reduction"
    exit 1
fi

echo "DEBUG: Prerequisites verified - matrix files found"

# ============================================================================
# R SCRIPT GENERATION FOR DIMENSIONALITY REDUCTION ANALYSIS
# ============================================================================
# Create comprehensive R script that performs the complete dimensionality
# reduction workflow including LSI, PCA, UMAP, and clustering analysis

cat > "$OUTPUT_DIR/dimensionality_reduction/${SAMPLE}_dim_reduction.R" << 'EOF'
#!/usr/bin/env Rscript
# ============================================================================
# DIMENSIONALITY REDUCTION PIPELINE FOR scATAC-seq DATA
# ============================================================================
# PURPOSE: Reduce dimensionality of sparse chromatin accessibility data
# METHODS: TF-IDF normalization → LSI → PCA → UMAP → K-means clustering
# OUTPUT: Low-dimensional embeddings and cluster assignments
# ============================================================================

# ============================================================================
# LIBRARY LOADING AND DEPENDENCY CHECK
# ============================================================================
# Load required packages for matrix operations, dimensionality reduction,
# clustering, and visualization. Each package serves specific functions:
# - Matrix: Sparse matrix operations for memory efficiency
# - irlba: Fast truncated SVD for LSI computation
# - uwot: UMAP implementation for non-linear embedding
# - cluster: Silhouette analysis for cluster quality assessment
# - ggplot2/dplyr: Data visualization and manipulation
# - RColorBrewer: Color palettes for cluster visualization
suppressPackageStartupMessages({
    library(Matrix)        # Sparse matrix operations
    library(irlba)         # Truncated SVD for LSI
    library(uwot)          # UMAP embedding
    library(Rtsne)         # t-SNE embedding (optional)
    library(cluster)       # Clustering and silhouette analysis
    library(ggplot2)       # Plotting and visualization
    library(dplyr)         # Data manipulation
    library(RColorBrewer)  # Color palettes
})

# ============================================================================
# COMMAND LINE ARGUMENT PROCESSING
# ============================================================================
# Handle sample name input from command line or environment variable
# This allows the script to be run independently or as part of the pipeline
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    # Fallback to environment variable or default sample
    sample_name <- Sys.getenv("SAMPLE", "R26-Nestin-Ctrl-adult")
} else {
    sample_name <- args[1]
}

cat("Processing sample:", sample_name, "\n")

# ============================================================================
# PATH CONFIGURATION AND FILE VALIDATION
# ============================================================================
# Set up input and output directory paths
# Input: 10X Genomics format files (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
# Output: Dimensionality reduction results and visualizations
matrix_dir <- "filtered_peak_bc_matrix_qc"
output_dir <- "dimensionality_reduction"

# ============================================================================
# DATA LOADING: 10X GENOMICS FORMAT PEAK-BARCODE MATRIX
# ============================================================================
# Load sparse matrix data in 10X Genomics format
# This format is memory-efficient for scATAC-seq data which is highly sparse
# (most peak-cell combinations have zero accessibility)
cat("Reading peak-barcode matrix...\n")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")

# Validate input file existence
# Critical check to prevent downstream errors
if(!file.exists(matrix_path)) {
    stop("Matrix file not found: ", matrix_path)
}

# Load the three components of 10X format:
# 1. matrix.mtx.gz: Sparse matrix in Matrix Market format (peak x cell)
# 2. features.tsv.gz: Peak coordinates (chr:start-end format)
# 3. barcodes.tsv.gz: Cell barcodes (unique identifiers per cell)
peak_matrix <- readMM(matrix_path)
features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Assign meaningful row and column names to the sparse matrix
# Rows: Peak coordinates (genomic regions)
# Columns: Cell barcodes (individual cells)
rownames(peak_matrix) <- features[,1]  # Peak coordinates (chr:start-end)
colnames(peak_matrix) <- barcodes[,1]  # Cell barcodes

cat("Matrix dimensions:", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# ============================================================================
# QUALITY CONTROL FILTERING
# ============================================================================
# Apply conservative filtering to remove low-quality cells and uninformative peaks
# This step is crucial for downstream dimensionality reduction quality
# Edge cases handled: cells with extreme peak counts, peaks present in very few cells
cat("Performing basic filtering...\n")

# Calculate accessibility statistics for quality assessment
# cell_peak_counts: Number of accessible peaks per cell (cell complexity)
# peak_cell_counts: Number of cells with accessible peak (peak prevalence)
cell_peak_counts <- Matrix::colSums(peak_matrix > 0)
peak_cell_counts <- Matrix::rowSums(peak_matrix > 0)

# ============================================================================
# DATA QUALITY ASSESSMENT AND DISTRIBUTION ANALYSIS
# ============================================================================
# Comprehensive statistics to understand data quality and inform filtering thresholds
# This analysis helps identify outlier cells and low-quality peaks
cat("Peak count distribution per cell:\n")
cat("  Min:", min(cell_peak_counts), "\n")
cat("  Max:", max(cell_peak_counts), "\n")
cat("  Median:", median(cell_peak_counts), "\n")
cat("  Mean:", mean(cell_peak_counts), "\n")
cat("  Quantiles (10%, 25%, 50%, 75%, 90%, 95%, 98%):\n")
print(quantile(cell_peak_counts, c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.98)))
cat("  Cells with >0 peaks:", sum(cell_peak_counts > 0), "\n")
cat("  Cells with >10 peaks:", sum(cell_peak_counts > 10), "\n")
cat("  Cells with >100 peaks:", sum(cell_peak_counts > 100), "\n")

# ============================================================================
# ADAPTIVE FILTERING THRESHOLD CALCULATION
# ============================================================================
# Calculate data-driven filtering thresholds to handle varying data quality
# Conservative approach to preserve biological signal while removing noise
# Edge case handling: Ensures minimum thresholds even for low-quality data

# Minimum peaks per cell: Focus on cells with meaningful accessibility profiles
# Use 100 as minimum threshold, but adapt if 95th percentile is higher
min_peaks_per_cell <- max(100, quantile(cell_peak_counts[cell_peak_counts > 0], 0.05))

# Maximum peaks per cell: 98th percentile to remove potential doublets (more stringent)
max_peaks_per_cell <- quantile(cell_peak_counts, 0.98)

# Minimum cells per peak: Remove very rare peaks that likely represent noise
# Use adaptive threshold: at least 0.1% of cells or minimum 10 cells
min_cells_per_peak <- max(10, ceiling(ncol(peak_matrix) * 0.001))

cat("Adjusted filtering thresholds:\n")
cat("  Min peaks per cell:", min_peaks_per_cell, "\n")
cat("  Max peaks per cell:", max_peaks_per_cell, "\n")
cat("  Min cells per peak:", min_cells_per_peak, "\n")

# ============================================================================
# FILTER APPLICATION AND MATRIX SUBSETTING
# ============================================================================
# Apply calculated thresholds to create filtered matrix for downstream analysis
# This removes low-quality cells and uninformative peaks

# Identify cells and peaks that pass quality filters
good_cells <- which(cell_peak_counts >= min_peaks_per_cell & 
                   cell_peak_counts <= max_peaks_per_cell)
good_peaks <- which(peak_cell_counts >= min_cells_per_peak)

# Create filtered matrix with only high-quality cells and informative peaks
filtered_matrix <- peak_matrix[good_peaks, good_cells]
cat("After filtering:", nrow(filtered_matrix), "peaks x", ncol(filtered_matrix), "cells\n")

# ============================================================================
# TF-IDF NORMALIZATION FOR scATAC-seq DATA
# ============================================================================
# Apply Term Frequency-Inverse Document Frequency normalization
# This is essential preprocessing for LSI in scATAC-seq data because:
# 1. Accounts for varying sequencing depth across cells (TF normalization)
# 2. Down-weights ubiquitous peaks while emphasizing cell-type-specific peaks (IDF)
# 3. Reduces technical noise and enhances biological signal
cat("Performing TF-IDF normalization...\n")

# ============================================================================
# MATRIX FORMAT CONVERSION FOR EFFICIENT COMPUTATION
# ============================================================================
# Ensure matrix is in compressed sparse column format for efficient operations
# This format is optimal for column-wise operations needed in TF-IDF calculation
if (class(filtered_matrix)[1] != "dgCMatrix") {
    filtered_matrix <- as(filtered_matrix, "CsparseMatrix")
}

# ============================================================================
# TERM FREQUENCY (TF) CALCULATION
# ============================================================================
# Normalize each cell by its total accessibility counts
# This accounts for differences in sequencing depth and cell size
# TF(peak, cell) = count(peak, cell) / total_counts(cell)
tf_matrix <- filtered_matrix
tf_matrix@x <- tf_matrix@x / rep.int(Matrix::colSums(tf_matrix), diff(tf_matrix@p))

# ============================================================================
# INVERSE DOCUMENT FREQUENCY (IDF) CALCULATION
# ============================================================================
# Calculate IDF to down-weight ubiquitous peaks and emphasize rare peaks
# IDF(peak) = log(total_cells / cells_with_peak)
# Higher IDF values indicate more cell-type-specific peaks
idf_values <- log(ncol(tf_matrix) / Matrix::rowSums(tf_matrix > 0))

# Handle edge case: peaks present in all cells (would cause infinite IDF)
idf_values[is.infinite(idf_values)] <- 0

# ============================================================================
# TF-IDF MATRIX CONSTRUCTION
# ============================================================================
# Apply IDF weighting to TF-normalized matrix
# Final TF-IDF(peak, cell) = TF(peak, cell) * IDF(peak)
# This creates a matrix where cell-type-specific accessible peaks have higher weights
tfidf_matrix <- tf_matrix
for(i in 1:nrow(tfidf_matrix)) {
    tfidf_matrix[i, ] <- tfidf_matrix[i, ] * idf_values[i]
}

# ============================================================================
# LATENT SEMANTIC INDEXING (LSI) - PRIMARY DIMENSIONALITY REDUCTION
# ============================================================================
# LSI is the gold standard for scATAC-seq dimensionality reduction because:
# 1. Captures major axes of chromatin accessibility variation
# 2. Reduces noise while preserving biological signal
# 3. Creates interpretable components representing cell states
# 4. Handles sparse, high-dimensional scATAC-seq data effectively
cat("Performing LSI...\n")

# ============================================================================
# COMPONENT NUMBER CALCULATION WITH SAFETY CHECKS
# ============================================================================
# Determine optimal number of LSI components while avoiding numerical issues
# Conservative approach to prevent SVD convergence problems
min_dim <- min(nrow(tfidf_matrix), ncol(tfidf_matrix))
max_safe_components <- min_dim - 3  # Conservative buffer for numerical stability

# ============================================================================
# ADAPTIVE COMPONENT SELECTION
# ============================================================================
# Balance between capturing biological variation and computational efficiency
# Edge case handling: Ensures valid component numbers for small datasets

# Start with reasonable default (50 components) but adapt to data size
n_components <- min(50, max_safe_components)
n_components <- max(2, n_components)  # Ensure minimum of 2 components

# Final safety validation to prevent SVD errors
if (n_components >= min_dim) {
    n_components <- max(2, min_dim - 3)
}

cat("Matrix dimensions:", nrow(tfidf_matrix), "x", ncol(tfidf_matrix), "\n")
cat("Min dimension:", min_dim, ", Max safe components:", max_safe_components, "\n")
cat("Using", n_components, "LSI components\n")

# ============================================================================
# PRE-SVD VALIDATION
# ============================================================================
# Critical validation to prevent runtime errors in SVD computation
if (n_components >= min_dim || n_components < 1) {
    stop("Invalid n_components: ", n_components, " for matrix of size ", 
         nrow(tfidf_matrix), "x", ncol(tfidf_matrix))
}

# ============================================================================
# SINGULAR VALUE DECOMPOSITION (SVD) FOR LSI
# ============================================================================
# Use irlba for efficient truncated SVD computation
# This is much faster than full SVD for large sparse matrices
# Set seed for reproducible results across runs
set.seed(42)
lsi_result <- irlba(tfidf_matrix, nv = n_components, nu = n_components)

# ============================================================================
# LSI EMBEDDING EXTRACTION AND FORMATTING
# ============================================================================
# Extract right singular vectors (V matrix) as cell embeddings
# Skip first component (LSI_1) as it often correlates with sequencing depth
# This is a standard practice in scATAC-seq analysis
lsi_embeddings <- lsi_result$v[, 2:n_components]
rownames(lsi_embeddings) <- colnames(filtered_matrix)  # Cell barcodes
colnames(lsi_embeddings) <- paste0("LSI_", 1:(n_components-1))  # Component names

cat("LSI completed. Using", ncol(lsi_embeddings), "components\n")

# ============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA) ON LSI COMPONENTS
# ============================================================================
# PCA on LSI components provides additional dimensionality reduction and:
# 1. Further reduces noise in the LSI space
# 2. Creates orthogonal components for downstream analysis
# 3. Enables variance explained calculations for component selection
# 4. Provides standardized coordinates for clustering and visualization
cat("Performing PCA on LSI components...\n")

# ============================================================================
# PCA COMPUTATION WITH CENTERING (NO SCALING)
# ============================================================================
# Center data but don't scale (LSI components already normalized)
# Centering removes mean shifts between components
# No scaling preserves relative importance of LSI components
pca_result <- prcomp(lsi_embeddings, center = TRUE, scale. = FALSE)

# ============================================================================
# PCA EMBEDDING EXTRACTION
# ============================================================================
# Extract all principal components for downstream analysis
# PCA embeddings maintain full dimensionality of LSI space
pca_embeddings <- pca_result$x

# ============================================================================
# VARIANCE EXPLAINED CALCULATION
# ============================================================================
# Calculate proportion of variance explained by each PC
# This helps determine how many components to retain for downstream analysis
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(variance_explained)

# ============================================================================
# PCA SUMMARY STATISTICS
# ============================================================================
# Report variance explained by first 10 components (or all if fewer)
# This provides insight into data complexity and component importance
cat("PCA completed. First 10 components explain", 
    round(cumulative_variance[min(10, length(cumulative_variance))] * 100, 2), 
    "% of variance\n")

# ============================================================================
# UNIFORM MANIFOLD APPROXIMATION AND PROJECTION (UMAP)
# ============================================================================
# UMAP creates 2D visualization of high-dimensional LSI space while:
# 1. Preserving local neighborhood structure (similar cells stay close)
# 2. Maintaining global topology (distant cell types remain separated)
# 3. Enabling intuitive visualization of cell populations
# 4. Providing coordinates for interactive plotting and exploration
cat("Computing UMAP embedding...\n")

# ============================================================================
# UMAP PARAMETER CONFIGURATION
# ============================================================================
# n_neighbors=30: Balance between local and global structure preservation
# min_dist=0.3: Controls how tightly points are packed (0.1-0.5 typical)
# metric="cosine": Appropriate for normalized LSI components
# verbose=TRUE: Monitor progress for large datasets
# Seed ensures reproducible embeddings across runs
set.seed(42)
umap_result <- umap(lsi_embeddings[, 1:min(30, ncol(lsi_embeddings))], 
                   n_neighbors = 30, 
                   min_dist = 0.3,
                   metric = "cosine",
                   verbose = TRUE)

# ============================================================================
# UMAP EMBEDDING EXTRACTION AND FORMATTING
# ============================================================================
# Extract 2D coordinates and create data frame for downstream analysis
# Maintains cell barcode mapping for cluster assignment and visualization
umap_df <- data.frame(
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2],
    Barcode = rownames(lsi_embeddings)
)

# ============================================================================
# t-SNE EMBEDDING (OPTIONAL - COMMENTED FOR PERFORMANCE)
# ============================================================================
# t-SNE provides alternative non-linear dimensionality reduction
# Commented out by default due to computational cost and UMAP superiority
# Can be enabled for comparative analysis or specific use cases
# t-SNE parameters: perplexity adapted to cell count, dims=2 for visualization
# cat("Computing t-SNE embedding...\n")
# set.seed(42)
# tsne_result <- Rtsne(lsi_embeddings[, 1:min(30, ncol(lsi_embeddings))],
#                     dims = 2, perplexity = min(30, floor(nrow(lsi_embeddings)/5)),
#                     verbose = TRUE, check_duplicates = FALSE)
# 
# tsne_df <- data.frame(
#     tSNE1 = tsne_result$Y[, 1],
#     tSNE2 = tsne_result$Y[, 2],
#     Barcode = rownames(lsi_embeddings)
# )

# ============================================================================
# K-MEANS CLUSTERING ON LSI COMPONENTS
# ============================================================================
# Clustering identifies distinct cell populations based on chromatin accessibility
# K-means is chosen for its speed and effectiveness on LSI-reduced data
# Clusters represent putative cell types or states
cat("Performing clustering...\n")

# ============================================================================
# CLUSTER NUMBER SPECIFICATION
# ============================================================================
# Set number of clusters based on expected cell types in the experiment
# This can be adjusted based on biological knowledge or exploratory analysis
n_clusters <- 8  # Adjust based on expected cell types

# ============================================================================
# K-MEANS CLUSTERING EXECUTION
# ============================================================================
# Use first 20 LSI components (or fewer if available)
# nstart=20: Multiple random initializations for robust results
# Reproducible results ensured by setting seed
set.seed(42)
cluster_result <- kmeans(lsi_embeddings[, 1:min(20, ncol(lsi_embeddings))], 
                        centers = n_clusters, 
                        nstart = 20)

# ============================================================================
# CLUSTER ASSIGNMENT TO UMAP DATA FRAME
# ============================================================================
# Add cluster assignments to UMAP coordinates for visualization
# Convert to factor for proper categorical handling in plots
umap_df$Cluster <- factor(cluster_result$cluster)

# ============================================================================
# CLUSTER QUALITY ASSESSMENT - SILHOUETTE ANALYSIS
# ============================================================================
# Silhouette score measures cluster cohesion and separation
# Values: -1 (poor) to +1 (excellent clustering)
# Helps validate cluster number choice and overall quality
sil_score <- silhouette(cluster_result$cluster, 
                       dist(lsi_embeddings[, 1:min(10, ncol(lsi_embeddings))]))
avg_sil_width <- mean(sil_score[, "sil_width"])
cat("Average silhouette width:", round(avg_sil_width, 3), "\n")

# ============================================================================
# VISUALIZATION GENERATION
# ============================================================================
# Create publication-quality plots for dimensionality reduction results
# Visualizations help interpret clustering and assess data quality
cat("Generating plots...\n")

# ============================================================================
# UMAP CLUSTER VISUALIZATION
# ============================================================================
# Create scatter plot of UMAP embedding colored by cluster assignment
# Point styling: small size and transparency for large datasets
# Color palette: qualitative colors for distinct cluster identification
p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    labs(title = paste("UMAP Clustering -", sample_name),
         subtitle = paste("n =", nrow(umap_df), "cells,", n_clusters, "clusters")) +
    theme_minimal() +
    theme(legend.position = "right")

# ============================================================================
# SAVE UMAP CLUSTER PLOT
# ============================================================================
# High-resolution PNG for publication and presentation use
# Dimensions optimized for readability and detail preservation
ggsave(file.path(output_dir, paste0(sample_name, "_umap_clusters.png")), 
       p1, width = 10, height = 8, dpi = 300)

# ============================================================================
# PCA VARIANCE EXPLAINED ANALYSIS
# ============================================================================
# Calculate variance explained by each principal component
# Focus on first 20 components for interpretability
# Both individual and cumulative variance for component selection guidance
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2))[1:20]
pca_var_df <- data.frame(
    PC = 1:length(var_explained),
    Variance_Explained = var_explained,
    Cumulative = cumsum(var_explained)
)

# ============================================================================
# PCA VARIANCE EXPLAINED VISUALIZATION
# ============================================================================
# Dual-axis plot showing individual and cumulative variance explained
# Bar chart: individual PC contribution
# Line plot: cumulative variance for component selection guidance
p2 <- ggplot(pca_var_df, aes(x = PC, y = Variance_Explained)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", size = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(title = paste("PCA Variance Explained -", sample_name),
         x = "Principal Component", 
         y = "Proportion of Variance Explained") +
    theme_minimal() +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance"))

# ============================================================================
# SAVE PCA VARIANCE PLOT
# ============================================================================
# Export variance explained plot for component selection and quality assessment
ggsave(file.path(output_dir, paste0(sample_name, "_pca_variance.png")), 
       p2, width = 10, height = 6, dpi = 300)

# ============================================================================
# DATA EXPORT AND PERSISTENCE
# ============================================================================
# Save all computed embeddings and results for downstream analysis
# TSV format ensures compatibility with various analysis tools
cat("Saving results...\n")

# ============================================================================
# LSI EMBEDDINGS EXPORT
# ============================================================================
# Save LSI components for downstream integration and analysis
# Row names: cell barcodes, Column names: LSI component identifiers
write.table(lsi_embeddings, 
           file.path(output_dir, paste0(sample_name, "_lsi_embeddings.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# ============================================================================
# PCA EMBEDDINGS EXPORT
# ============================================================================
# Save PCA coordinates for alternative analysis approaches
# Maintains full dimensionality of PCA space
write.table(pca_embeddings, 
           file.path(output_dir, paste0(sample_name, "_pca_embeddings.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# ============================================================================
# UMAP EMBEDDINGS AND CLUSTERS EXPORT
# ============================================================================
# Save 2D UMAP coordinates with cluster assignments
# Essential for visualization and interactive analysis
write.table(umap_df, 
           file.path(output_dir, paste0(sample_name, "_umap_embeddings.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save cluster assignments
cluster_df <- data.frame(
    Barcode = names(cluster_result$cluster),
    Cluster = cluster_result$cluster,
    Silhouette_Width = sil_score[, "sil_width"]
)

write.table(cluster_df, 
           file.path(output_dir, paste0(sample_name, "_clusters.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save summary statistics
summary_stats <- data.frame(
    Metric = c("Original_peaks", "Original_cells", "Filtered_peaks", "Filtered_cells",
               "LSI_components", "PCA_components", "Clusters", "Avg_silhouette_width"),
    Value = c(nrow(peak_matrix), ncol(peak_matrix), nrow(filtered_matrix), 
              ncol(filtered_matrix), ncol(lsi_embeddings), ncol(pca_embeddings),
              n_clusters, avg_sil_width)
)

write.table(summary_stats, 
           file.path(output_dir, paste0(sample_name, "_dim_reduction_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save filtered matrix (optional, for downstream analysis)
Matrix::writeMM(filtered_matrix, 
               file.path(output_dir, paste0(sample_name, "_filtered_matrix.mtx")))

cat("Analysis complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("Files generated:\n")
cat("  - LSI embeddings:", paste0(sample_name, "_lsi_embeddings.tsv"), "\n")
cat("  - PCA embeddings:", paste0(sample_name, "_pca_embeddings.tsv"), "\n") 
cat("  - UMAP embeddings:", paste0(sample_name, "_umap_embeddings.tsv"), "\n")
cat("  - Cluster assignments:", paste0(sample_name, "_clusters.tsv"), "\n")
cat("  - Summary statistics:", paste0(sample_name, "_dim_reduction_summary.tsv"), "\n")
cat("  - Plots: UMAP clusters and PCA variance\n")
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/dimensionality_reduction/${SAMPLE}_dim_reduction.R"

# Run the R script
echo "DEBUG: Running dimensionality reduction analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    echo "DEBUG: You can install with: conda install -c conda-forge r-base"
    exit 1
fi

# Try to run the analysis
echo "DEBUG: Starting R analysis for $SAMPLE..."
Rscript "dimensionality_reduction/${SAMPLE}_dim_reduction.R" "$SAMPLE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: R analysis completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="dimensionality_reduction/${SAMPLE}_dim_reduction_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "Analysis Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: R analysis failed. Check the R script and dependencies."
    echo "Required R packages: Matrix, irlba, uwot, cluster, ggplot2, dplyr, RColorBrewer"
    echo "Install with: install.packages(c('Matrix', 'irlba', 'uwot', 'cluster', 'ggplot2', 'dplyr', 'RColorBrewer'))"
fi

echo "Output files expected in: $OUTPUT_DIR/dimensionality_reduction/"
echo "  - ${SAMPLE}_lsi_embeddings.tsv"
echo "  - ${SAMPLE}_pca_embeddings.tsv"
echo "  - ${SAMPLE}_umap_embeddings.tsv"
echo "  - ${SAMPLE}_clusters.tsv"
echo "  - ${SAMPLE}_dim_reduction_summary.tsv"
echo "  - ${SAMPLE}_umap_clusters.png"
echo "  - ${SAMPLE}_pca_variance.png"

echo "========================================="
echo "Step 9 complete for $SAMPLE"
echo "Dimensionality reduction analysis finished"
echo "End time: $(date)"
echo "========================================="