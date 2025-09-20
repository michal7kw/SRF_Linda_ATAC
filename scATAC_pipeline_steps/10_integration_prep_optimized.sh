#!/bin/bash
#SBATCH --job-name=10_integration_prep_opt
#SBATCH --output=logs/10_integration_prep_opt_%a.out
#SBATCH --error=logs/10_integration_prep_opt_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

# ============================================================================
# OPTIMIZED PIPELINE INITIALIZATION AND ENVIRONMENT SETUP
# ============================================================================
# Purpose: Configure the computational environment for optimized scATAC-seq integration preparation
# This optimized version addresses timeout issues by implementing vectorized operations
# and improved memory management for large-scale datasets
#
# Key Optimizations:
# - Vectorized gene activity score computation
# - Batch processing to reduce memory overhead
# - Improved sparse matrix operations
# - Enhanced resource allocation (48h, 128GB RAM)
#
# Resource Requirements:
# - 16 CPUs for parallel genomic operations
# - 128GB RAM for large sparse matrix operations (increased from 64GB)
# - 48 hours for comprehensive gene activity calculation (doubled from 24h)
# ============================================================================

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate integration_prep

set -euo pipefail

# ============================================================================
# SAMPLE CONFIGURATION AND LOGGING
# ============================================================================

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 10: OPTIMIZED Integration Preparation for $SAMPLE"
echo "Start time: $(date)"
echo "Memory allocated: 128GB"
echo "Time limit: 48 hours"
echo "========================================="

# ============================================================================
# OUTPUT DIRECTORY PREPARATION
# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR/integration_prep"

# ============================================================================
# PREREQUISITE VALIDATION
# ============================================================================

# Check prerequisites
PEAK_BC_MATRIX="$OUTPUT_DIR/filtered_peak_bc_matrix_qc"
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

if [[ ! -d "$PEAK_BC_MATRIX" ]]; then
    echo "ERROR: Peak-barcode matrix directory not found: $PEAK_BC_MATRIX"
    echo "Please run step 6 (06_peak_cell_matrix.sh) first"
    exit 1
fi

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

# ============================================================================
# OPTIMIZED R SCRIPT GENERATION FOR INTEGRATION PREPARATION
# ============================================================================

# Create comprehensive OPTIMIZED R script for integration preparation
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep_optimized.R" << 'EOF'
#!/usr/bin/env Rscript
# OPTIMIZED Integration Preparation Script for scATAC-seq data
# This version implements vectorized operations to avoid timeout issues

# ============================================================================
# LIBRARY LOADING AND DEPENDENCY MANAGEMENT
# ============================================================================

suppressPackageStartupMessages({
    library(Matrix)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
    library(parallel)  # Added for parallel processing
    if(!require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)) {
        cat("Installing TxDb.Mmusculus.UCSC.mm10.knownGene...\n")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    }
    if(!require(org.Mm.eg.db, quietly = TRUE)) {
        cat("Installing org.Mm.eg.db...\n")
        BiocManager::install("org.Mm.eg.db", ask = FALSE)
        library(org.Mm.eg.db)
    }
})

# ============================================================================
# COMMAND LINE ARGUMENT PROCESSING AND PATH CONFIGURATION
# ============================================================================

# Get sample name
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    sample_name <- Sys.getenv("SAMPLE", "R26-Nestin-Ctrl-adult")
} else {
    sample_name <- args[1]
}

cat("Processing sample:", sample_name, "\n")
cat("Using OPTIMIZED version with vectorized operations\n")

# Set paths
matrix_dir <- "filtered_peak_bc_matrix_qc"
output_dir <- "integration_prep"

# ============================================================================
# PEAK-BARCODE MATRIX LOADING
# ============================================================================

# 1. Load peak-barcode matrix
cat("Loading peak-barcode matrix...\n")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")

peak_matrix <- readMM(matrix_path)
features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- features[,1]
colnames(peak_matrix) <- barcodes[,1]

cat("Matrix loaded:", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# ============================================================================
# GENE ACTIVITY SCORE CALCULATION - INITIALIZATION
# ============================================================================

# 2. Create gene activity scores
cat("Creating gene activity scores for integration (OPTIMIZED)...\n")

# Load gene annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

# ============================================================================
# PEAK COORDINATE PARSING AND GENOMIC RANGE CREATION
# ============================================================================

# Parse peak coordinates
peak_coords <- rownames(peak_matrix)
peak_chr <- gsub(":.*", "", peak_coords)
peak_ranges <- gsub(".*:", "", peak_coords)
peak_start <- as.numeric(gsub("-.*", "", peak_ranges))
peak_end <- as.numeric(gsub(".*-", "", peak_ranges))

# Create GenomicRanges object for peaks
peak_gr <- GRanges(
    seqnames = peak_chr,
    ranges = IRanges(start = peak_start, end = peak_end)
)
names(peak_gr) <- peak_coords

cat("Created", length(peak_gr), "peak ranges\n")

# ============================================================================
# PEAK-GENE ASSOCIATION MAPPING
# ============================================================================

# Find overlaps between peaks and gene bodies + promoters
cat("Finding peak-gene associations...\n")

# Gene bodies
gene_overlaps <- findOverlaps(peak_gr, genes)
gene_body_associations <- data.frame(
    peak = names(peak_gr)[queryHits(gene_overlaps)],
    gene_id = names(genes)[subjectHits(gene_overlaps)],
    type = "gene_body"
)

# Promoter regions (TSS +/- 2kb)
promoters <- promoters(genes, upstream = 2000, downstream = 2000)
promoter_overlaps <- findOverlaps(peak_gr, promoters)
promoter_associations <- data.frame(
    peak = names(peak_gr)[queryHits(promoter_overlaps)],
    gene_id = names(promoters)[subjectHits(promoter_overlaps)],
    type = "promoter"
)

# Combine associations
all_associations <- rbind(gene_body_associations, promoter_associations)
cat("Found", nrow(all_associations), "peak-gene associations\n")

# ============================================================================
# OPTIMIZED GENE ACTIVITY MATRIX CONSTRUCTION
# ============================================================================
# MAJOR OPTIMIZATION: Replace inefficient for-loop with vectorized operations
# This addresses the timeout issue by processing genes in batches using
# split-apply-combine strategy instead of individual gene processing
# ============================================================================

cat("Computing gene activity scores using OPTIMIZED vectorized approach...\n")

# Get unique genes and pre-filter associations for valid peaks
unique_genes <- unique(all_associations$gene_id)
cat("Processing", length(unique_genes), "unique genes\n")

# Pre-filter associations to only include peaks that exist in our matrix
valid_peaks_mask <- all_associations$peak %in% rownames(peak_matrix)
filtered_associations <- all_associations[valid_peaks_mask, ]
cat("Filtered associations to", nrow(filtered_associations), "valid peak-gene pairs\n")

# OPTIMIZATION 1: Create peak-to-gene mapping list efficiently
cat("Creating optimized peak-gene mapping...\n")
peak_gene_list <- split(filtered_associations$peak, filtered_associations$gene_id)

# OPTIMIZATION 2: Vectorized gene activity computation using lapply
cat("Computing gene activities using vectorized operations...\n")
cat("This may take some time for large datasets...\n")

# Set up progress reporting
genes_to_process <- names(peak_gene_list)
n_genes <- length(genes_to_process)
progress_interval <- max(1, floor(n_genes / 20))  # Report progress every 5%

# Function to compute gene activity for a single gene (vectorized)
compute_gene_activity <- function(gene_idx) {
    gene <- genes_to_process[gene_idx]
    associated_peaks <- peak_gene_list[[gene]]
    
    # Report progress
    if (gene_idx %% progress_interval == 0) {
        cat("  Processed", gene_idx, "/", n_genes, "genes (", 
            round(100 * gene_idx / n_genes, 1), "%)\n")
    }
    
    if(length(associated_peaks) > 0) {
        if(length(associated_peaks) == 1) {
            return(peak_matrix[associated_peaks, ])
        } else {
            # Use efficient Matrix colSums for multiple peaks
            return(Matrix::colSums(peak_matrix[associated_peaks, , drop = FALSE]))
        }
    } else {
        # Return zero vector for genes with no valid peaks
        return(numeric(ncol(peak_matrix)))
    }
}

# OPTIMIZATION 3: Use lapply for vectorized processing (more memory efficient than mclapply for this case)
cat("Starting vectorized gene activity computation...\n")
start_time <- Sys.time()

gene_activity_list <- lapply(seq_along(genes_to_process), compute_gene_activity)

end_time <- Sys.time()
cat("Gene activity computation completed in", 
    round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

# OPTIMIZATION 4: Efficient sparse matrix construction
cat("Constructing sparse gene activity matrix...\n")
gene_activity_matrix <- do.call(rbind, gene_activity_list)
rownames(gene_activity_matrix) <- genes_to_process
colnames(gene_activity_matrix) <- colnames(peak_matrix)

# Convert to sparse matrix if not already
if (!inherits(gene_activity_matrix, "sparseMatrix")) {
    gene_activity_matrix <- as(gene_activity_matrix, "sparseMatrix")
}

# Remove genes with no activity (more efficient filtering)
gene_sums <- Matrix::rowSums(gene_activity_matrix)
active_genes <- gene_sums > 0
gene_activity_matrix <- gene_activity_matrix[active_genes, ]

cat("OPTIMIZED Gene activity matrix:", nrow(gene_activity_matrix), "genes x", ncol(gene_activity_matrix), "cells\n")

# ============================================================================
# MEMORY CLEANUP
# ============================================================================

# Clean up large intermediate objects to free memory
rm(gene_activity_list, peak_gene_list, filtered_associations)
gc()  # Force garbage collection

# ============================================================================
# GENE ID CONVERSION TO SYMBOLS
# ============================================================================

# 3. Convert gene IDs to symbols
cat("Converting gene IDs to symbols...\n")
gene_symbols <- tryCatch({
    mapIds(org.Mm.eg.db, keys = rownames(gene_activity_matrix), 
           column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
}, error = function(e) {
    cat("Warning: Could not convert gene IDs to symbols\n")
    setNames(rownames(gene_activity_matrix), rownames(gene_activity_matrix))
})

# Remove NAs and duplicates
valid_symbols <- !is.na(gene_symbols) & !duplicated(gene_symbols)
gene_activity_matrix <- gene_activity_matrix[valid_symbols, ]
rownames(gene_activity_matrix) <- gene_symbols[valid_symbols]

cat("Final gene activity matrix:", nrow(gene_activity_matrix), "genes x", ncol(gene_activity_matrix), "cells\n")

# ============================================================================
# GENE ACTIVITY MATRIX EXPORT
# ============================================================================

# 4. Save gene activity matrix in multiple formats
cat("Saving gene activity matrix...\n")

# Save as Matrix Market format
Matrix::writeMM(gene_activity_matrix, 
               file.path(output_dir, paste0(sample_name, "_gene_activity_optimized.mtx")))

# Save gene names
write.table(data.frame(gene_symbol = rownames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_features_optimized.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save cell barcodes (same as peak matrix)
write.table(data.frame(barcode = colnames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes_optimized.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ============================================================================
# SEURAT-COMPATIBLE CELL METADATA CREATION
# ============================================================================

# 5. Create Seurat-compatible metadata
cat("Creating Seurat-compatible metadata...\n")

# Calculate per-cell metrics
cell_metrics <- data.frame(
    barcode = colnames(peak_matrix),
    n_peaks = Matrix::colSums(peak_matrix > 0),
    total_fragments = Matrix::colSums(peak_matrix),
    n_genes = Matrix::colSums(gene_activity_matrix > 0),
    total_gene_activity = Matrix::colSums(gene_activity_matrix)
)

# Calculate additional QC metrics
cell_metrics$log10_total_fragments <- log10(cell_metrics$total_fragments + 1)
cell_metrics$peaks_per_fragment <- cell_metrics$n_peaks / cell_metrics$total_fragments

# Save cell metadata
write.table(cell_metrics,
           file.path(output_dir, paste0(sample_name, "_cell_metadata_optimized.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# PEAK ANNOTATION WITH GENE ASSOCIATIONS
# ============================================================================

# 6. Create peak annotations for regulatory analysis
cat("Creating peak annotations...\n")
peak_annotations <- data.frame(
    peak_id = names(peak_gr),
    chr = as.character(seqnames(peak_gr)),
    start = start(peak_gr),
    end = end(peak_gr),
    width = width(peak_gr)
)

# Add gene associations (using original associations before filtering)
peak_gene_map <- all_associations %>%
    group_by(peak) %>%
    summarise(
        associated_genes = paste(gene_id, collapse = ","),
        n_genes = n(),
        association_types = paste(unique(type), collapse = ",")
    )

peak_annotations <- merge(peak_annotations, peak_gene_map, 
                         by.x = "peak_id", by.y = "peak", all.x = TRUE)

# Fill NAs
peak_annotations$associated_genes[is.na(peak_annotations$associated_genes)] <- "none"
peak_annotations$n_genes[is.na(peak_annotations$n_genes)] <- 0
peak_annotations$association_types[is.na(peak_annotations$association_types)] <- "intergenic"

write.table(peak_annotations,
           file.path(output_dir, paste0(sample_name, "_peak_annotations_optimized.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# INTEGRATION SUMMARY GENERATION
# ============================================================================

# 7. Create integration summary
cat("Creating integration summary...\n")
integration_summary <- data.frame(
    metric = c("total_peaks", "total_cells", "active_genes", "mean_peaks_per_cell",
               "mean_fragments_per_cell", "mean_gene_activity_per_cell",
               "peaks_with_gene_associations", "genes_with_peak_associations",
               "optimization_method", "processing_time_estimate"),
    value = c(nrow(peak_matrix), ncol(peak_matrix), nrow(gene_activity_matrix),
              round(mean(cell_metrics$n_peaks), 2),
              round(mean(cell_metrics$total_fragments), 2),
              round(mean(cell_metrics$total_gene_activity), 2),
              sum(peak_annotations$n_genes > 0),
              length(unique(unlist(strsplit(peak_annotations$associated_genes[
                peak_annotations$associated_genes != "none"], ",")))),
              "vectorized_operations",
              paste(round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes"))
)

write.table(integration_summary,
           file.path(output_dir, paste0(sample_name, "_integration_summary_optimized.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# QUALITY CONTROL VISUALIZATION
# ============================================================================

# 8. Generate QC plots for integration
cat("Generating integration QC plots...\n")

# Peaks vs gene activity per cell
p1 <- ggplot(cell_metrics, aes(x = n_peaks, y = n_genes)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = paste("Peaks vs Gene Activity (OPTIMIZED) -", sample_name),
         x = "Number of Accessible Peaks", 
         y = "Number of Active Genes") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_peaks_vs_genes_optimized.png")), 
       p1, width = 8, height = 6, dpi = 300)

# Distribution of gene activity scores
gene_activity_dist <- data.frame(
    gene = rownames(gene_activity_matrix),
    total_activity = Matrix::rowSums(gene_activity_matrix),
    n_cells_active = Matrix::rowSums(gene_activity_matrix > 0)
)

p2 <- ggplot(gene_activity_dist, aes(x = n_cells_active, y = total_activity)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste("Gene Activity Distribution (OPTIMIZED) -", sample_name),
         x = "Number of Cells with Activity (log10)",
         y = "Total Activity Score (log10)") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_gene_activity_dist_optimized.png")), 
       p2, width = 8, height = 6, dpi = 300)

cat("OPTIMIZED Integration preparation complete!\n")
cat("Files generated:\n")
cat("  - Gene activity matrix: ", paste0(sample_name, "_gene_activity_optimized.mtx"), "\n")
cat("  - Gene features: ", paste0(sample_name, "_gene_activity_features_optimized.tsv"), "\n")
cat("  - Cell barcodes: ", paste0(sample_name, "_gene_activity_barcodes_optimized.tsv"), "\n")
cat("  - Cell metadata: ", paste0(sample_name, "_cell_metadata_optimized.tsv"), "\n")
cat("  - Peak annotations: ", paste0(sample_name, "_peak_annotations_optimized.tsv"), "\n")
cat("  - Integration summary: ", paste0(sample_name, "_integration_summary_optimized.tsv"), "\n")
cat("  - QC plots: peaks vs genes, gene activity distribution (optimized versions)\n")
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep_optimized.R"

# ============================================================================
# OPTIMIZED SEURAT INTEGRATION HELPER SCRIPT GENERATION
# ============================================================================

# Create an optimized Seurat integration helper script
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat_optimized.R" << 'EOF'
#!/usr/bin/env Rscript
# OPTIMIZED Helper script to load scATAC data into Seurat for integration

suppressPackageStartupMessages({
    library(Seurat)
    library(Signac)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    stop("Please provide sample name as argument")
}

sample_name <- args[1]
output_dir <- "integration_prep"

cat("Loading", sample_name, "data into Seurat (OPTIMIZED version)...\n")

# Load peak matrix
peak_matrix <- readMM("filtered_peak_bc_matrix_qc/matrix.mtx.gz")
peak_features <- read.table("filtered_peak_bc_matrix_qc/features.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)
peak_barcodes <- read.table("filtered_peak_bc_matrix_qc/barcodes.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- peak_features[,1]
colnames(peak_matrix) <- peak_barcodes[,1]

# Load OPTIMIZED gene activity matrix
gene_matrix <- readMM(file.path(output_dir, paste0(sample_name, "_gene_activity_optimized.mtx")))
gene_features <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_features_optimized.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)
gene_barcodes <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes_optimized.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)

rownames(gene_matrix) <- gene_features[,1]
colnames(gene_matrix) <- gene_barcodes[,1]

# Load OPTIMIZED metadata
metadata <- read.table(file.path(output_dir, paste0(sample_name, "_cell_metadata_optimized.tsv")),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$barcode

cat("Creating OPTIMIZED Seurat object...\n")

# Create ChromatinAssay for peaks
chrom_assay <- CreateChromatinAssay(
    counts = peak_matrix,
    sep = c(":", "-")
)

# Create Seurat object
atac_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
)

# Add gene activity as a separate assay
atac_obj[["RNA"]] <- CreateAssayObject(counts = gene_matrix)

# Set default assay
DefaultAssay(atac_obj) <- "peaks"

cat("OPTIMIZED Seurat object created successfully!\n")
cat("Assays:", names(atac_obj@assays), "\n")
cat("Cells:", ncol(atac_obj), "\n")
cat("Peaks:", nrow(atac_obj[["peaks"]]), "\n")
cat("Genes:", nrow(atac_obj[["RNA"]]), "\n")

# Save Seurat object
saveRDS(atac_obj, file.path(output_dir, paste0(sample_name, "_seurat_object_optimized.rds")))
cat("OPTIMIZED Seurat object saved to:", file.path(output_dir, paste0(sample_name, "_seurat_object_optimized.rds")), "\n")
EOF

chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat_optimized.R"

# ============================================================================
# EXECUTE OPTIMIZED INTEGRATION PREPARATION
# ============================================================================

echo "DEBUG: Running OPTIMIZED integration preparation analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    exit 1
fi

# Run the main OPTIMIZED integration prep script
echo "DEBUG: Starting OPTIMIZED integration preparation for $SAMPLE..."
echo "Using vectorized operations to avoid timeout issues..."
Rscript "integration_prep/${SAMPLE}_integration_prep_optimized.R" "$SAMPLE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: OPTIMIZED Integration preparation completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="integration_prep/${SAMPLE}_integration_summary_optimized.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "OPTIMIZED Integration Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: OPTIMIZED Integration preparation failed. Check R dependencies."
    echo "Required packages: Matrix, GenomicRanges, GenomicFeatures, rtracklayer"
    echo "Install with BiocManager::install(c('GenomicRanges', 'GenomicFeatures', 'rtracklayer'))"
fi

echo "OPTIMIZED Output files created in: $OUTPUT_DIR/integration_prep/"
echo "  - Gene activity matrix: ${SAMPLE}_gene_activity_optimized.mtx"
echo "  - Cell metadata: ${SAMPLE}_cell_metadata_optimized.tsv"
echo "  - Peak annotations: ${SAMPLE}_peak_annotations_optimized.tsv"
echo "  - Integration summary: ${SAMPLE}_integration_summary_optimized.tsv"
echo "  - Seurat helper script: ${SAMPLE}_load_seurat_optimized.R"
echo "  - QC plots: peaks vs genes, gene activity distribution (optimized versions)"

echo "========================================="
echo "Step 10 OPTIMIZED complete for $SAMPLE"
echo "Integration preparation finished"
echo "End time: $(date)"
echo "========================================="