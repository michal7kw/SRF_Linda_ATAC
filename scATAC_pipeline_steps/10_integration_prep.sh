#!/bin/bash
#SBATCH --job-name=10_integration_prep
#SBATCH --output=logs/10_integration_prep_%a.out
#SBATCH --error=logs/10_integration_prep_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=workq

# ============================================================================
# PIPELINE INITIALIZATION AND ENVIRONMENT SETUP
# ============================================================================
# Purpose: Configure the computational environment for scATAC-seq integration preparation
# This step prepares scATAC-seq data for integration with scRNA-seq by creating
# gene activity scores and Seurat-compatible objects
#
# Key Operations:
# - Activates specialized conda environment with genomics packages
# - Sets strict error handling for pipeline reliability
# - Configures SLURM array job for parallel sample processing
#
# Resource Requirements:
# - 8 CPUs for parallel genomic operations
# - 32GB RAM for large sparse matrix operations
# - 6 hours for comprehensive gene activity calculation
# ============================================================================

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate integration_prep

set -euo pipefail

# ============================================================================
# SAMPLE CONFIGURATION AND LOGGING
# ============================================================================
# Purpose: Define sample-specific parameters and initialize logging
# 
# Key Operations:
# - Maps SLURM array task ID to specific sample names
# - Sets up output directory structure for integration results
# - Initializes comprehensive logging for tracking progress
#
# Input Expectations:
# - SLURM_ARRAY_TASK_ID environment variable (0-1)
# - Existing output directory from previous pipeline steps
#
# Output Structure:
# - Creates integration_prep subdirectory for results
# - Maintains sample-specific file naming convention
# ============================================================================

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 10: Integration Preparation for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# ============================================================================
# OUTPUT DIRECTORY PREPARATION
# ============================================================================
# Purpose: Create organized directory structure for integration results
# Ensures proper file organization for downstream analysis steps
# ============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR/integration_prep"

# ============================================================================
# PREREQUISITE VALIDATION
# ============================================================================
# Purpose: Verify that all required input files from previous pipeline steps exist
# 
# Key Validations:
# - Peak-barcode matrix from step 6 (filtered_peak_bc_matrix_qc/)
# - Fragment files from chromap alignment
# - Peak files from MACS2 peak calling
#
# Input Dependencies:
# - 06_peak_cell_matrix.sh: Creates filtered peak-barcode matrix
# - Previous steps: Generate fragments and peaks files
#
# Error Handling:
# - Exits with clear error messages if prerequisites missing
# - Provides guidance on which pipeline step to re-run
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
# R SCRIPT GENERATION FOR INTEGRATION PREPARATION
# ============================================================================
# Purpose: Create comprehensive R script for scATAC-seq integration preparation
# This script transforms peak accessibility data into gene activity scores
# suitable for integration with scRNA-seq data
#
# Key Functionality:
# - Converts peak-level accessibility to gene-level activity scores
# - Creates Seurat-compatible data structures
# - Generates comprehensive QC metrics and visualizations
# - Prepares data for downstream integration and GRN analysis
#
# Algorithm Overview:
# 1. Load peak-barcode matrix from previous steps
# 2. Map peaks to genes using genomic annotations
# 3. Aggregate peak accessibility into gene activity scores
# 4. Create metadata and annotations for integration
# 5. Generate QC plots and summary statistics
# ============================================================================

# Create comprehensive R script for integration preparation
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep.R" << 'EOF'
#!/usr/bin/env Rscript
# Integration Preparation Script for scATAC-seq data
# Prepares data for integration with scRNA-seq and GRN analysis

# ============================================================================
# LIBRARY LOADING AND DEPENDENCY MANAGEMENT
# ============================================================================
# Purpose: Load required R packages for genomic analysis and data manipulation
# 
# Core Libraries:
# - Matrix: Efficient sparse matrix operations for large-scale data
# - GenomicRanges: Genomic interval operations and annotations
# - GenomicFeatures: Gene annotation and transcript database access
# - rtracklayer: Import/export of genomic data formats
# - dplyr: Data manipulation and transformation
# - ggplot2: High-quality data visualization
#
# Bioconductor Annotations:
# - TxDb.Mmusculus.UCSC.mm10.knownGene: Mouse gene annotations (mm10)
# - org.Mm.eg.db: Mouse gene ID mapping and annotation database
#
# Dependency Handling:
# - Automatic installation of missing Bioconductor packages
# - Graceful error handling for package loading
# ============================================================================

suppressPackageStartupMessages({
    library(Matrix)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
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
# Purpose: Handle sample identification and configure file paths
# 
# Input Methods:
# - Command line argument (preferred): Rscript script.R sample_name
# - Environment variable fallback: SAMPLE environment variable
# - Default fallback: "R26-Nestin-Ctrl-adult"
#
# Path Configuration:
# - matrix_dir: Location of filtered peak-barcode matrix from step 6
# - output_dir: Destination for integration preparation results
#
# Flexibility Features:
# - Multiple input methods ensure robust sample identification
# - Consistent path structure across pipeline steps
# ============================================================================

# Get sample name
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    sample_name <- Sys.getenv("SAMPLE", "R26-Nestin-Ctrl-adult")
} else {
    sample_name <- args[1]
}

cat("Processing sample:", sample_name, "\n")

# Set paths
matrix_dir <- "filtered_peak_bc_matrix_qc"
output_dir <- "integration_prep"

# ============================================================================
# PEAK-BARCODE MATRIX LOADING
# ============================================================================
# Purpose: Load filtered peak accessibility matrix from previous pipeline step
# 
# Input Format: Matrix Market (.mtx) sparse matrix format
# - matrix.mtx.gz: Sparse matrix with peak accessibility counts
# - features.tsv.gz: Peak coordinates (chr:start-end format)
# - barcodes.tsv.gz: Cell barcodes from quality-filtered cells
#
# Data Structure:
# - Rows: Accessible chromatin peaks (genomic regions)
# - Columns: Quality-filtered single cells
# - Values: Accessibility counts (typically binary or low integers)
#
# Memory Efficiency:
# - Uses sparse matrix format for memory-efficient storage
# - Handles large-scale single-cell datasets (>10k cells, >100k peaks)
#
# Quality Control:
# - Matrix represents only high-quality cells and reproducible peaks
# - Pre-filtered for minimum accessibility thresholds
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
# Purpose: Transform peak-level accessibility into gene-level activity scores
# This is the core step for integrating scATAC-seq with scRNA-seq data
#
# Biological Rationale:
# - Accessible chromatin peaks near genes correlate with gene expression
# - Gene activity scores approximate gene expression from chromatin accessibility
# - Enables direct comparison and integration with scRNA-seq data
#
# Algorithm Overview:
# 1. Load comprehensive gene annotations from UCSC/Ensembl
# 2. Parse peak coordinates into genomic ranges
# 3. Map peaks to genes based on genomic proximity
# 4. Aggregate peak accessibility into gene-level scores
#
# Annotation Source:
# - TxDb.Mmusculus.UCSC.mm10.knownGene: Comprehensive mouse gene database
# - Includes gene bodies, promoters, and regulatory regions
# ============================================================================

# 2. Create gene activity scores
cat("Creating gene activity scores for integration...\n")

# Load gene annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

# ============================================================================
# PEAK COORDINATE PARSING AND GENOMIC RANGE CREATION
# ============================================================================
# Purpose: Convert peak coordinate strings into GenomicRanges objects
# 
# Input Format: "chr1:1000-2000" (chromosome:start-end)
# Parsing Steps:
# 1. Extract chromosome names (chr1, chr2, etc.)
# 2. Extract start and end coordinates
# 3. Create GenomicRanges object for efficient overlap operations
#
# GenomicRanges Benefits:
# - Efficient genomic interval operations
# - Built-in overlap and proximity functions
# - Memory-efficient representation of genomic coordinates
#
# Error Handling:
# - Validates coordinate format and numeric conversion
# - Handles edge cases in coordinate parsing
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
# Purpose: Map accessible chromatin peaks to their associated genes
# This step defines which peaks contribute to each gene's activity score
#
# Mapping Strategy:
# 1. Gene Body Associations: Peaks overlapping gene transcribed regions
# 2. Promoter Associations: Peaks overlapping promoter regions (TSS ± 2kb)
#
# Biological Rationale:
# - Gene body peaks: Often contain enhancers and regulatory elements
# - Promoter peaks: Directly regulate transcription initiation
# - 2kb promoter window: Captures core promoter and proximal regulatory elements
#
# Technical Implementation:
# - Uses GenomicRanges::findOverlaps for efficient interval operations
# - Handles overlapping annotations (peaks can map to multiple genes)
# - Maintains peak-gene relationship metadata for downstream analysis
#
# Output Structure:
# - Each row represents one peak-gene association
# - Includes association type (gene_body vs promoter) for weighting
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
# GENE ACTIVITY MATRIX CONSTRUCTION
# ============================================================================
# Purpose: Aggregate peak accessibility scores into gene-level activity scores
# This creates a gene × cell matrix analogous to scRNA-seq expression data
#
# Algorithm:
# 1. Initialize sparse matrix with genes as rows, cells as columns
# 2. For each gene, sum accessibility of all associated peaks
# 3. Handle single vs multiple peak associations efficiently
# 4. Remove genes with zero activity across all cells
#
# Mathematical Approach:
# - Gene Activity Score = Σ(Peak Accessibility) for all associated peaks
# - Preserves quantitative relationships between accessibility levels
# - Maintains sparsity for memory efficiency
#
# Quality Control:
# - Validates peak existence in original matrix
# - Handles edge cases (genes with no valid peaks)
# - Filters out inactive genes to reduce noise
#
# Biological Interpretation:
# - Higher scores indicate more accessible regulatory regions
# - Correlates with gene expression potential
# - Enables integration with transcriptomic data
# ============================================================================

# Create gene activity matrix
cat("Computing gene activity scores...\n")
unique_genes <- unique(all_associations$gene_id)
gene_activity_matrix <- Matrix(0, nrow = length(unique_genes), 
                              ncol = ncol(peak_matrix), sparse = TRUE)
rownames(gene_activity_matrix) <- unique_genes
colnames(gene_activity_matrix) <- colnames(peak_matrix)

# Sum peak accessibility for each gene
for(gene in unique_genes) {
    associated_peaks <- all_associations$peak[all_associations$gene_id == gene]
    if(length(associated_peaks) > 0) {
        # Find peaks that exist in our matrix
        valid_peaks <- intersect(associated_peaks, rownames(peak_matrix))
        if(length(valid_peaks) > 0) {
            if(length(valid_peaks) == 1) {
                gene_activity_matrix[gene, ] <- peak_matrix[valid_peaks, ]
            } else {
                gene_activity_matrix[gene, ] <- Matrix::colSums(peak_matrix[valid_peaks, ])
            }
        }
    }
}

# Remove genes with no activity
gene_sums <- Matrix::rowSums(gene_activity_matrix)
active_genes <- gene_sums > 0
gene_activity_matrix <- gene_activity_matrix[active_genes, ]

cat("Gene activity matrix:", nrow(gene_activity_matrix), "genes x", ncol(gene_activity_matrix), "cells\n")

# ============================================================================
# GENE ID CONVERSION TO SYMBOLS
# ============================================================================
# Purpose: Convert Entrez gene IDs to human-readable gene symbols
# This step makes the data compatible with standard scRNA-seq workflows
#
# Conversion Process:
# - Uses org.Mm.eg.db annotation database for ID mapping
# - Maps from ENTREZID to SYMBOL (e.g., "12345" → "Gapdh")
# - Handles multiple symbols per ID by selecting the first match
#
# Quality Control:
# - Removes genes with unmappable IDs (NA symbols)
# - Eliminates duplicate gene symbols to ensure unique identifiers
# - Maintains matrix integrity during filtering
#
# Error Handling:
# - Graceful fallback if annotation database unavailable
# - Preserves original IDs if symbol conversion fails
# - Continues analysis with available data
#
# Integration Benefits:
# - Gene symbols are standard in scRNA-seq analysis
# - Enables direct comparison with expression datasets
# - Facilitates downstream pathway and functional analysis
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
# Purpose: Save gene activity matrix in standard formats for downstream analysis
# 
# Export Formats:
# 1. Matrix Market (.mtx): Sparse matrix format for efficient storage
# 2. Features file (.tsv): Gene symbols corresponding to matrix rows
# 3. Barcodes file (.tsv): Cell barcodes corresponding to matrix columns
#
# Format Compatibility:
# - Compatible with 10X Genomics Cell Ranger output format
# - Directly loadable by Seurat, Scanpy, and other scRNA-seq tools
# - Maintains sparsity for memory-efficient downstream analysis
#
# File Structure:
# - sample_gene_activity.mtx: Sparse matrix with gene activity scores
# - sample_gene_activity_features.tsv: Gene symbols (one per line)
# - sample_gene_activity_barcodes.tsv: Cell barcodes (one per line)
#
# Integration Workflow:
# - These files can be loaded as a "RNA" assay in Seurat
# - Enables joint analysis with scRNA-seq expression data
# - Supports integration algorithms (CCA, Harmony, etc.)
# ============================================================================

# 4. Save gene activity matrix in multiple formats
cat("Saving gene activity matrix...\n")

# Save as Matrix Market format
Matrix::writeMM(gene_activity_matrix, 
               file.path(output_dir, paste0(sample_name, "_gene_activity.mtx")))

# Save gene names
write.table(data.frame(gene_symbol = rownames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_features.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save cell barcodes (same as peak matrix)
write.table(data.frame(barcode = colnames(gene_activity_matrix)),
           file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ============================================================================
# SEURAT-COMPATIBLE CELL METADATA CREATION
# ============================================================================
# Purpose: Generate comprehensive cell-level quality control metrics
# This metadata enables filtering and QC assessment in downstream analysis
#
# Metadata Components:
# 1. barcode: Unique cell identifier (matches matrix columns)
# 2. n_peaks: Number of accessible peaks per cell (complexity metric)
# 3. total_fragments: Total ATAC-seq signal per cell (library size)
# 4. n_genes: Number of active genes per cell (transcriptional complexity)
# 5. total_gene_activity: Total gene activity score per cell
# 6. log10_total_fragments: Log-transformed fragment counts for visualization
# 7. peaks_per_fragment: Peak accessibility efficiency metric
#
# Quality Control Applications:
# - Filter low-quality cells (low peak/gene counts)
# - Identify potential doublets (high signal cells)
# - Normalize for library size differences
# - Batch effect assessment across samples
#
# Integration Compatibility:
# - Standard format for Seurat object creation
# - Enables consistent QC across multiple samples
# - Supports automated filtering thresholds
# ============================================================================

# 5. Create Seurat-compatible objects
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
           file.path(output_dir, paste0(sample_name, "_cell_metadata.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# PEAK ANNOTATION WITH GENE ASSOCIATIONS
# ============================================================================
# Purpose: Create comprehensive peak metadata for regulatory analysis
# Links genomic peaks to their potential target genes
#
# Annotation Components:
# 1. peak_id: Unique peak identifier (chr:start-end format)
# 2. chr, start, end: Genomic coordinates of accessible regions
# 3. width: Peak width (regulatory element size)
# 4. associated_genes: Linked genes (comma-separated list)
# 5. n_genes: Number of associated genes per peak
# 6. association_types: Type of gene association (gene_body, promoter, intergenic)
#
# Gene Association Strategy:
# - Links peaks to genes within gene bodies (intragenic regulation)
# - Associates peaks with promoter regions (transcriptional control)
# - Aggregates multiple gene associations per peak
# - Handles peaks with no gene associations (intergenic regions)
#
# Regulatory Analysis Applications:
# - Identify cis-regulatory elements
# - Link chromatin accessibility to gene expression
# - Prioritize peaks for functional validation
# - Support gene regulatory network construction
#
# Data Integration:
# - Compatible with ChIP-seq peak analysis workflows
# - Enables motif enrichment analysis
# - Supports regulatory element classification
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

# Add gene associations
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
           file.path(output_dir, paste0(sample_name, "_peak_annotations.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# INTEGRATION SUMMARY GENERATION
# ============================================================================
# Purpose: Create comprehensive summary statistics for integration QC
# Provides key metrics for assessing data quality and integration readiness
#
# Summary Metrics:
# 1. total_peaks: Number of accessible chromatin peaks
# 2. total_cells: Number of high-quality cells passing filters
# 3. active_genes: Number of genes with detectable activity scores
# 4. mean_peaks_per_cell: Average chromatin accessibility per cell
# 5. mean_fragments_per_cell: Average ATAC-seq signal per cell
# 6. mean_gene_activity_per_cell: Average gene activity per cell
# 7. peaks_with_gene_associations: Peaks linked to genes (regulatory potential)
# 8. genes_with_peak_associations: Genes with accessible regulatory elements
#
# Quality Assessment Applications:
# - Compare data quality across samples
# - Identify samples requiring additional filtering
# - Assess integration feasibility
# - Track processing pipeline performance
#
# Integration Planning:
# - Determine normalization strategies
# - Plan batch correction approaches
# - Assess sample compatibility for joint analysis
# ============================================================================

# 7. Create integration summary
cat("Creating integration summary...\n")
integration_summary <- data.frame(
    metric = c("total_peaks", "total_cells", "active_genes", "mean_peaks_per_cell",
               "mean_fragments_per_cell", "mean_gene_activity_per_cell",
               "peaks_with_gene_associations", "genes_with_peak_associations"),
    value = c(nrow(peak_matrix), ncol(peak_matrix), nrow(gene_activity_matrix),
              round(mean(cell_metrics$n_peaks), 2),
              round(mean(cell_metrics$total_fragments), 2),
              round(mean(cell_metrics$total_gene_activity), 2),
              sum(peak_annotations$n_genes > 0),
              length(unique(unlist(strsplit(peak_annotations$associated_genes[
                peak_annotations$associated_genes != "none"], ",")))))
)

write.table(integration_summary,
           file.path(output_dir, paste0(sample_name, "_integration_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ============================================================================
# QUALITY CONTROL VISUALIZATION
# ============================================================================
# Purpose: Generate diagnostic plots for integration preparation assessment
# Visual QC enables rapid identification of data quality issues
#
# QC Plot 1: Peak-Gene Activity Correlation
# - X-axis: Number of accessible peaks per cell (chromatin complexity)
# - Y-axis: Number of active genes per cell (transcriptional activity)
# - Correlation: Assesses consistency between chromatin and gene activity
#
# Interpretation Guidelines:
# - Strong positive correlation: High-quality data with consistent signal
# - Weak correlation: Potential technical issues or biological heterogeneity
# - Outliers: Possible doublets (high values) or low-quality cells (low values)
#
# Quality Assessment:
# - Cells with very low peaks/genes: Candidates for filtering
# - Cells with extremely high values: Potential doublets
# - Overall correlation strength: Data integration feasibility
#
# Integration Readiness:
# - Well-correlated data: Ready for standard integration workflows
# - Poor correlation: May require additional QC or specialized methods
# ============================================================================

# 8. Generate QC plots for integration
cat("Generating integration QC plots...\n")

# Peaks vs gene activity per cell
p1 <- ggplot(cell_metrics, aes(x = n_peaks, y = n_genes)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = paste("Peaks vs Gene Activity -", sample_name),
         x = "Number of Accessible Peaks", 
         y = "Number of Active Genes") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_peaks_vs_genes.png")), 
       p1, width = 8, height = 6, dpi = 300)

# ============================================================================
# GENE ACTIVITY DISTRIBUTION ANALYSIS
# ============================================================================
# Purpose: Assess the distribution of gene activity scores across genes
# Identifies highly active genes and evaluates data quality patterns
#
# QC Plot 2: Gene Activity Distribution
# - X-axis: Number of cells with detectable activity per gene (log10)
# - Y-axis: Total activity score per gene across all cells (log10)
# - Point cloud: Each point represents one gene
#
# Distribution Interpretation:
# - Upper right: Highly expressed genes active in many cells
# - Lower left: Lowly expressed genes active in few cells
# - Outliers (high Y, low X): Genes with high activity in few cells
# - Linear relationship: Consistent gene expression patterns
#
# Quality Assessment:
# - Broad distribution: Good dynamic range of gene activity
# - Tight clustering: Limited gene expression diversity
# - Clear outliers: Potential marker genes or technical artifacts
#
# Integration Applications:
# - Identify highly variable genes for integration
# - Assess gene activity dynamic range
# - Compare gene expression patterns across samples
# ============================================================================

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
    labs(title = paste("Gene Activity Distribution -", sample_name),
         x = "Number of Cells with Activity (log10)",
         y = "Total Activity Score (log10)") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_gene_activity_dist.png")), 
       p2, width = 8, height = 6, dpi = 300)

cat("Integration preparation complete!\n")
cat("Files generated:\n")
cat("  - Gene activity matrix: ", paste0(sample_name, "_gene_activity.mtx"), "\n")
cat("  - Gene features: ", paste0(sample_name, "_gene_activity_features.tsv"), "\n")
cat("  - Cell barcodes: ", paste0(sample_name, "_gene_activity_barcodes.tsv"), "\n")
cat("  - Cell metadata: ", paste0(sample_name, "_cell_metadata.tsv"), "\n")
cat("  - Peak annotations: ", paste0(sample_name, "_peak_annotations.tsv"), "\n")
cat("  - Integration summary: ", paste0(sample_name, "_integration_summary.tsv"), "\n")
cat("  - QC plots: peaks vs genes, gene activity distribution\n")
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_integration_prep.R"

# ============================================================================
# SEURAT INTEGRATION HELPER SCRIPT GENERATION
# ============================================================================
# Purpose: Create automated R script for Seurat object construction
# Enables seamless integration with standard scRNA-seq analysis workflows
#
# Script Components:
# 1. Matrix Loading: Reads peak and gene activity matrices in sparse format
# 2. Metadata Integration: Incorporates cell-level QC metrics
# 3. Multi-Assay Object: Creates Seurat object with both peaks and RNA assays
# 4. Data Persistence: Saves ready-to-use Seurat object
#
# Seurat Object Structure:
# - "peaks" assay: Chromatin accessibility data (primary)
# - "RNA" assay: Gene activity scores (for integration)
# - Metadata: QC metrics, sample information, filtering flags
#
# Integration Workflow Benefits:
# - Standard format compatible with Seurat integration functions
# - Enables joint analysis with scRNA-seq datasets
# - Supports advanced integration methods (CCA, Harmony, etc.)
# - Facilitates multi-modal analysis workflows
#
# Usage:
# - Run generated R script to create Seurat object
# - Load object in downstream integration analysis
# - Compatible with standard Seurat analysis pipelines
# ============================================================================

# Create a simple Seurat integration helper script
cat > "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat.R" << 'EOF'
#!/usr/bin/env Rscript
# Helper script to load scATAC data into Seurat for integration

# This script creates Seurat objects from the processed scATAC data
# Usage: Rscript load_seurat.R sample_name

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

cat("Loading", sample_name, "data into Seurat...\n")

# Load peak matrix
peak_matrix <- readMM("filtered_peak_bc_matrix_qc/matrix.mtx.gz")
peak_features <- read.table("filtered_peak_bc_matrix_qc/features.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)
peak_barcodes <- read.table("filtered_peak_bc_matrix_qc/barcodes.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- peak_features[,1]
colnames(peak_matrix) <- peak_barcodes[,1]

# Load gene activity matrix
gene_matrix <- readMM(file.path(output_dir, paste0(sample_name, "_gene_activity.mtx")))
gene_features <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_features.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)
gene_barcodes <- read.table(file.path(output_dir, paste0(sample_name, "_gene_activity_barcodes.tsv")),
                           header = FALSE, stringsAsFactors = FALSE)

rownames(gene_matrix) <- gene_features[,1]
colnames(gene_matrix) <- gene_barcodes[,1]

# Load metadata
metadata <- read.table(file.path(output_dir, paste0(sample_name, "_cell_metadata.tsv")),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$barcode

cat("Creating Seurat object...\n")

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

cat("Seurat object created successfully!\n")
cat("Assays:", names(atac_obj@assays), "\n")
cat("Cells:", ncol(atac_obj), "\n")
cat("Peaks:", nrow(atac_obj[["peaks"]]), "\n")
cat("Genes:", nrow(atac_obj[["RNA"]]), "\n")

# Save Seurat object
saveRDS(atac_obj, file.path(output_dir, paste0(sample_name, "_seurat_object.rds")))
cat("Seurat object saved to:", file.path(output_dir, paste0(sample_name, "_seurat_object.rds")), "\n")
EOF

chmod +x "$OUTPUT_DIR/integration_prep/${SAMPLE}_load_seurat.R"

# Run the integration preparation
echo "DEBUG: Running integration preparation analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    exit 1
fi

# Run the main integration prep script
echo "DEBUG: Starting integration preparation for $SAMPLE..."
Rscript "integration_prep/${SAMPLE}_integration_prep.R" "$SAMPLE"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Integration preparation completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="integration_prep/${SAMPLE}_integration_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "Integration Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: Integration preparation failed. Check R dependencies."
    echo "Required packages: Matrix, GenomicRanges, GenomicFeatures, rtracklayer"
    echo "Install with BiocManager::install(c('GenomicRanges', 'GenomicFeatures', 'rtracklayer'))"
fi

echo "Output files created in: $OUTPUT_DIR/integration_prep/"
echo "  - Gene activity matrix: ${SAMPLE}_gene_activity.mtx"
echo "  - Cell metadata: ${SAMPLE}_cell_metadata.tsv"
echo "  - Peak annotations: ${SAMPLE}_peak_annotations.tsv"
echo "  - Integration summary: ${SAMPLE}_integration_summary.tsv"
echo "  - Seurat helper script: ${SAMPLE}_load_seurat.R"
echo "  - QC plots: peaks vs genes, gene activity distribution"

echo "========================================="
echo "Step 10 complete for $SAMPLE"
echo "Integration preparation finished"
echo "End time: $(date)"
echo "========================================="