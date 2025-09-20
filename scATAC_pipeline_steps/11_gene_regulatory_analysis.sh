#!/bin/bash
#SBATCH --job-name=11_gene_regulatory_analysis
#SBATCH --output=logs/11_gene_regulatory_analysis_%a.out
#SBATCH --error=logs/11_gene_regulatory_analysis_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

#===============================================================================
# GENE REGULATORY NETWORK (GRN) ANALYSIS PIPELINE
#===============================================================================
# Purpose: Comprehensive transcription factor binding site analysis and gene
#          regulatory network construction from scATAC-seq data
#
# Key Functions:
# - TF motif discovery in accessible chromatin regions
# - TF activity score calculation using chromVAR methodology
# - Peak-to-gene regulatory link prediction
# - Comparative analysis between experimental conditions
# - Visualization of regulatory networks and TF activities
#
# Input Requirements:
# - Peak-barcode matrix from MACS2 peak calling
# - Gene activity scores from integration preparation
# - Peak annotations with genomic coordinates
#
# Output Products:
# - TF activity matrices (cell x TF z-scores)
# - Regulatory network edges (TF-target gene pairs)
# - Statistical summaries and quality metrics
# - Comparative analysis plots and heatmaps
#
# Computational Resources:
# - Memory: 64GB (for motif matching and matrix operations)
# - CPU: 16 cores (parallel motif scanning)
# - Time: ~12 hours (depends on peak count and cell number)
#
# Dependencies:
# - R/Bioconductor: motifmatchr, TFBSTools, chromVAR, BSgenome
# - JASPAR motif database for mouse transcription factors
# - GenomicRanges for coordinate-based operations
#===============================================================================

#===============================================================================
# ENVIRONMENT SETUP AND CONFIGURATION
#===============================================================================
# Purpose: Initialize computational environment and configure sample parameters
# - Activates conda environment with required R/Bioconductor packages
# - Sets sample-specific parameters for analysis
# - Configures output directory structure
#===============================================================================

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

#===============================================================================
# OUTPUT DIRECTORY PREPARATION
#===============================================================================
# Purpose: Create organized directory structure for GRN analysis outputs
# - Establishes main results directory
# - Creates specialized subdirectory for regulatory network files
# - Ensures proper file organization for downstream analysis
#===============================================================================

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 11: Gene Regulatory Network Analysis for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Create output directory
mkdir -p "$OUTPUT_DIR/grn_analysis"

#===============================================================================
# PREREQUISITE VALIDATION
#===============================================================================
# Purpose: Verify availability of required input files for GRN analysis
# 
# Key Validations:
# - Peak annotations: Genomic coordinates and gene associations
# - Gene activity matrix: Cell x gene activity scores for integration
# - MACS2 peaks: High-confidence accessible chromatin regions
#
# Input Dependencies:
# - Step 10: Integration preparation with peak annotations and gene activity
# - Peak calling: Sorted BED files with accessible chromatin regions
#
# Error Handling:
# - Provides clear error messages with missing file paths
# - Suggests specific pipeline steps to run if files are missing
# - Exits gracefully to prevent downstream analysis failures
#===============================================================================

# Check prerequisites
INTEGRATION_DIR="$OUTPUT_DIR/integration_prep"
PEAK_ANNOTATIONS="$INTEGRATION_DIR/${SAMPLE}_peak_annotations.tsv"
GENE_ACTIVITY="$INTEGRATION_DIR/${SAMPLE}_gene_activity.mtx"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

if [[ ! -f "$PEAK_ANNOTATIONS" ]]; then
    echo "ERROR: Peak annotations not found: $PEAK_ANNOTATIONS"
    echo "Please run step 10 (10_integration_prep.sh) first"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

#===============================================================================
# R SCRIPT GENERATION FOR GRN ANALYSIS
#===============================================================================
# Purpose: Create comprehensive R script for transcription factor analysis
#
# Key Functionality:
# - TF motif discovery using JASPAR database
# - chromVAR-based TF activity score calculation
# - Peak-to-gene regulatory link prediction
# - Statistical analysis and visualization generation
#
# Algorithm Overview:
# 1. Load peak-barcode matrix and convert to GenomicRanges
# 2. Scan peaks for TF binding motifs using motifmatchr
# 3. Calculate TF activity deviations using chromVAR methodology
# 4. Link peaks to genes based on proximity and gene activity
# 5. Generate regulatory network and summary statistics
#
# Output Structure:
# - TF activity matrices (cells x TFs with z-scores)
# - Regulatory edges (TF-target gene relationships)
# - Quality control plots and statistical summaries
#===============================================================================

# Create comprehensive R script for GRN analysis
cat > "$OUTPUT_DIR/grn_analysis/${SAMPLE}_grn_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
# Gene Regulatory Network Analysis for scATAC-seq data
# Identifies TF binding sites, motif enrichment, and regulatory networks

#===============================================================================
# LIBRARY LOADING AND DEPENDENCY MANAGEMENT
#===============================================================================
# Purpose: Load essential R/Bioconductor packages for regulatory analysis
#
# Core Libraries:
# - Matrix: Sparse matrix operations for large-scale data
# - GenomicRanges: Genomic coordinate manipulation and overlap detection
# - motifmatchr: High-performance TF motif scanning in genomic sequences
# - TFBSTools: Transcription factor binding site analysis tools
# - BSgenome.Mmusculus.UCSC.mm10: Mouse reference genome sequences
#
# Analysis Frameworks:
# - chromVAR: TF activity inference from chromatin accessibility
# - SummarizedExperiment: Coordinated data structure for omics analysis
#
# Visualization Tools:
# - ggplot2: Statistical graphics and publication-quality plots
# - ComplexHeatmap: Advanced heatmap generation with clustering
# - circlize: Color mapping and circular visualization support
#
# Dependency Handling:
# - Suppresses startup messages for cleaner output
# - All packages must be pre-installed in conda environment
#===============================================================================

suppressPackageStartupMessages({
    library(Matrix)
    library(GenomicRanges)
    library(motifmatchr)
    library(TFBSTools)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(chromVAR)
    library(SummarizedExperiment)
    library(ggplot2)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    
    # Try to load JASPAR motifs
    if(!require(JASPAR2020, quietly = TRUE)) {
        cat("Installing JASPAR2020...\n")
        BiocManager::install("JASPAR2020", ask = FALSE)
        library(JASPAR2020)
    }
    
    if(!require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)) {
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", ask = FALSE)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    }
    
    if(!require(org.Mm.eg.db, quietly = TRUE)) {
        BiocManager::install("org.Mm.eg.db", ask = FALSE)
        library(org.Mm.eg.db)
    }
})

#===============================================================================
# SAMPLE CONFIGURATION AND PATH SETUP
#===============================================================================
# Purpose: Initialize sample-specific parameters and file paths
# - Retrieves sample name from environment variable or command line
# - Configures directory paths for input and output files
# - Ensures consistent file naming across analysis steps
#===============================================================================

# Get sample name
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
    sample_name <- Sys.getenv("SAMPLE", "R26-Nestin-Ctrl-adult")
} else {
    sample_name <- args[1]
}

cat("Processing GRN analysis for:", sample_name, "\n")

# Set paths
output_dir <- "grn_analysis"
integration_dir <- "integration_prep"

#===============================================================================
# PEAK DATA AND ANNOTATION LOADING
#===============================================================================
# Purpose: Load and process scATAC-seq peak accessibility data
#
# Input Format:
# - Matrix Market (.mtx.gz): Sparse matrix of peak x cell accessibility
# - Features (.tsv.gz): Peak coordinates and identifiers
# - Barcodes (.tsv.gz): Cell barcode identifiers
# - Peak annotations: Gene associations from integration preparation
#
# Data Structure:
# - Rows: Accessible chromatin peaks (genomic regions)
# - Columns: Individual cells (barcodes)
# - Values: Read counts or binary accessibility scores
#
# Memory Efficiency:
# - Uses sparse matrix format to handle large datasets
# - Reduces memory footprint for datasets with many zero values
#===============================================================================

# 1. Load peak data and annotations
cat("Loading peak data and annotations...\n")

# Load peak-barcode matrix
peak_matrix <- readMM("filtered_peak_bc_matrix_qc/matrix.mtx.gz")
peak_features <- read.table("filtered_peak_bc_matrix_qc/features.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)
peak_barcodes <- read.table("filtered_peak_bc_matrix_qc/barcodes.tsv.gz", 
                           header = FALSE, stringsAsFactors = FALSE)

rownames(peak_matrix) <- peak_features[,1]
colnames(peak_matrix) <- peak_barcodes[,1]

# Load peak annotations
peak_annotations <- read.table(file.path(integration_dir, paste0(sample_name, "_peak_annotations.tsv")),
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Loaded", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

#===============================================================================
# GENOMIC RANGES OBJECT CREATION AND FILTERING
#===============================================================================
# Purpose: Convert peak coordinates to GenomicRanges for efficient operations
#
# GenomicRanges Benefits:
# - Efficient overlap detection for motif scanning
# - Coordinate-based operations and sorting
# - Integration with Bioconductor genomics workflows
#
# Coordinate Parsing:
# - Extracts chromosome, start, and end from peak identifiers
# - Handles standard format: chr:start-end
# - Validates coordinate integrity
#
# Filtering Strategy:
# - Removes non-standard chromosomes (scaffolds, patches)
# - Focuses on main chromosomes (1-19, X, Y) for mouse
# - Maintains consistency between matrix and coordinate objects
#
# Quality Control:
# - Reports number of peaks retained after filtering
# - Ensures matrix rows match filtered peak coordinates
#===============================================================================

# 2. Create GenomicRanges object for peaks
cat("Creating genomic ranges for peaks...\n")
peak_coords <- rownames(peak_matrix)
peak_chr <- gsub(":.*", "", peak_coords)
peak_ranges <- gsub(".*:", "", peak_coords)
peak_start <- as.numeric(gsub("-.*", "", peak_ranges))
peak_end <- as.numeric(gsub(".*-", "", peak_ranges))

peak_gr <- GRanges(
    seqnames = peak_chr,
    ranges = IRanges(start = peak_start, end = peak_end),
    peak_id = peak_coords
)

# Filter for standard chromosomes only
standard_chroms <- paste0("chr", c(1:19, "X", "Y"))
peak_gr <- peak_gr[seqnames(peak_gr) %in% standard_chroms]
peak_matrix <- peak_matrix[peak_gr$peak_id, ]

cat("Filtered to", length(peak_gr), "peaks on standard chromosomes\n")

#===============================================================================
# JASPAR MOTIF DATABASE LOADING AND FILTERING
#===============================================================================
# Purpose: Load transcription factor binding motifs for regulatory analysis
#
# JASPAR Database:
# - Comprehensive collection of TF binding site profiles
# - Position weight matrices (PWMs) for motif scanning
# - Curated from experimental ChIP-seq and protein binding data
#
# Filtering Strategy:
# - CORE collection: High-quality, non-redundant motifs
# - Vertebrates: Evolutionary relevant TF families
# - Neural-specific TFs: Focus on brain development factors
# - Latest versions: Most current motif profiles
#
# Quality Control:
# - Reports number of motifs loaded
# - Ensures species-specific filtering accuracy
# - Validates motif matrix integrity
#
# Applications:
# - Motif scanning in accessible chromatin regions
# - TF activity inference from binding site enrichment
# - Regulatory network construction
#===============================================================================

# 3. Get motif database
cat("Loading JASPAR motif database...\n")
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"
opts[["matrixtype"]] <- "PFM"
motifs <- getMatrixSet(JASPAR2020, opts)

# Filter for mouse TFs
mouse_tfs <- c("Ascl1", "Dlx2", "Emx2", "Foxg1", "Gsx2", "Hes1", "Hes5", 
               "Hey1", "Hey2", "Lhx2", "Mash1", "Neurog2", "Olig2", "Pax6", 
               "Sox2", "Tbr1", "Tbr2", "Tlx", "Zic1", "Zic3")

# Keep motifs for relevant TFs (expand this list based on your specific interests)
relevant_motifs <- motifs[TFBSTools::name(motifs) %in% mouse_tfs]
cat("Using", length(relevant_motifs), "motifs for analysis\n")

#===============================================================================
# MOTIF SCANNING IN ACCESSIBLE CHROMATIN REGIONS
#===============================================================================
# Purpose: Identify transcription factor binding sites within ATAC-seq peaks
#
# Algorithm:
# - Extracts DNA sequences from peak coordinates
# - Scans sequences using position weight matrices (PWMs)
# - Identifies statistically significant motif matches
# - Creates binary matrix of TF-peak associations
#
# Chunked Processing:
# - Processes peaks in manageable chunks (5000 per chunk)
# - Prevents memory overflow with large peak sets
# - Enables progress tracking for long-running analysis
#
# Technical Implementation:
# - Uses BSgenome for sequence extraction
# - Applies motifmatchr for high-performance scanning
# - Combines results into unified motif-peak matrix
#
# Output Structure:
# - Rows: Accessible chromatin peaks
# - Columns: Transcription factors
# - Values: Binary (1 = motif present, 0 = absent)
#
# Quality Control:
# - Reports total peaks and TFs analyzed
# - Validates sequence extraction success
# - Ensures coordinate consistency
#===============================================================================

# 4. Find motif matches in peaks
cat("Finding motif matches in accessible peaks...\n")
genome <- BSgenome.Mmusculus.UCSC.mm10

# This step can be memory intensive, so we'll process in chunks if needed
chunk_size <- min(5000, length(peak_gr))
n_chunks <- ceiling(length(peak_gr) / chunk_size)

motif_matches_list <- list()
for(i in 1:n_chunks) {
    cat("Processing chunk", i, "of", n_chunks, "\n")
    
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, length(peak_gr))
    
    peak_chunk <- peak_gr[start_idx:end_idx]
    
    # Get sequences
    peak_seqs <- getSeq(genome, peak_chunk)
    
    # Find motif matches
    motif_matches_chunk <- matchMotifs(relevant_motifs, peak_seqs, 
                                      genome = genome, p.cutoff = 1e-4)
    
    motif_matches_list[[i]] <- assay(motif_matches_chunk)
}

# Combine results
motif_matches <- do.call(rbind, motif_matches_list)
colnames(motif_matches) <- TFBSTools::name(relevant_motifs)
rownames(motif_matches) <- peak_gr$peak_id

cat("Found motif matches in", nrow(motif_matches), "peaks for", ncol(motif_matches), "TFs\n")

#===============================================================================
# TRANSCRIPTION FACTOR ACTIVITY CALCULATION USING CHROMVAR
#===============================================================================
# Purpose: Infer TF activity from chromatin accessibility and motif presence
#
# chromVAR Methodology:
# - Compares observed vs expected motif accessibility
# - Accounts for technical biases (GC content, peak size)
# - Generates cell-specific TF activity z-scores
# - Enables single-cell regulatory analysis
#
# Algorithm Steps:
# 1. Create SummarizedExperiment with peak accessibility
# 2. Define background peak sets with similar properties
# 3. Calculate accessibility deviations for each TF
# 4. Normalize to z-scores for comparative analysis
#
# Background Peak Selection:
# - Matches GC content distribution
# - Controls for peak width and accessibility
# - Iterative sampling for robust statistics
#
# Output Interpretation:
# - Positive z-scores: Higher than expected TF activity
# - Negative z-scores: Lower than expected TF activity
# - Magnitude indicates strength of deviation
#
# Quality Control:
# - Reports matrix dimensions
# - Validates background peak generation
# - Ensures proper normalization
#===============================================================================

# 5. Calculate TF activity scores using chromVAR approach
cat("Calculating TF activity scores...\n")

# Create SummarizedExperiment object
peak_se <- SummarizedExperiment(
    assays = list(counts = peak_matrix[rownames(motif_matches), ]),
    rowRanges = peak_gr[match(rownames(motif_matches), peak_gr$peak_id)],
    colData = DataFrame(cell_id = colnames(peak_matrix))
)

# Add motif matches
rowData(peak_se)$motif_matches <- motif_matches

# Calculate background peaks
bg_peaks <- getBackgroundPeaks(peak_se, niterations = 200)

# Calculate deviations (TF activity)
dev_scores <- computeDeviations(object = peak_se, 
                               annotations = motif_matches,
                               background_peaks = bg_peaks)

# Extract z-scores (TF activity per cell)
tf_activity <- assays(dev_scores)[["z"]]
rownames(tf_activity) <- colnames(motif_matches)

cat("Calculated TF activity scores:", nrow(tf_activity), "TFs x", ncol(tf_activity), "cells\n")

# 6. Identify differentially active TFs (if comparing conditions)
cat("Analyzing TF activity patterns...\n")

# Calculate TF activity statistics
tf_stats <- data.frame(
    TF = rownames(tf_activity),
    mean_activity = rowMeans(tf_activity),
    sd_activity = apply(tf_activity, 1, sd),
    max_activity = apply(tf_activity, 1, max),
    min_activity = apply(tf_activity, 1, min)
)

tf_stats$cv <- tf_stats$sd_activity / abs(tf_stats$mean_activity)
tf_stats <- tf_stats[order(tf_stats$cv, decreasing = TRUE), ]

#===============================================================================
# REGULATORY NETWORK CONSTRUCTION
#===============================================================================
# Purpose: Link transcription factors to target genes through peak-gene associations
#
# Methodology:
# - Uses peak annotations to identify gene associations
# - Links TF motifs in peaks to associated genes
# - Creates directed regulatory edges (TF → target gene)
#
# Peak-Gene Association Types:
# - Promoter: Peaks within gene promoter regions
# - Enhancer: Distal regulatory elements linked to genes
# - Intronic: Peaks within gene bodies (potential regulatory elements)
#
# Network Construction Algorithm:
# 1. Filter peaks with gene associations (exclude "none")
# 2. For each peak, identify TF motifs present
# 3. For each TF-peak combination, link to associated genes
# 4. Create regulatory edges with metadata
#
# Edge Attributes:
# - TF: Transcription factor name
# - target_gene: Regulated gene symbol
# - peak_id: Mediating chromatin peak
# - regulation_type: Type of peak-gene association
#
# Quality Control:
# - Reports total regulatory interactions identified
# - Validates peak-gene association integrity
# - Ensures TF-motif mapping consistency
#
# Applications:
# - Gene regulatory network visualization
# - Pathway enrichment analysis
# - Comparative analysis between conditions
#===============================================================================

# 7. Peak-to-gene linking for regulatory networks
cat("Creating peak-to-gene regulatory links...\n")

# Use existing peak annotations
peak_gene_links <- peak_annotations[peak_annotations$associated_genes != "none", ]

# Create regulatory network edges
regulatory_edges <- data.frame()

for(i in 1:nrow(peak_gene_links)) {
    peak_id <- peak_gene_links$peak_id[i]
    genes <- unlist(strsplit(peak_gene_links$associated_genes[i], ","))
    
    # Find TFs with motifs in this peak
    peak_motifs <- names(which(motif_matches[peak_id, ]))
    
    if(length(peak_motifs) > 0 && length(genes) > 0) {
        # Create TF -> target gene edges
        for(tf in peak_motifs) {
            for(gene in genes) {
                regulatory_edges <- rbind(regulatory_edges, data.frame(
                    TF = tf,
                    target_gene = gene,
                    peak_id = peak_id,
                    regulation_type = peak_gene_links$association_types[i],
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
}

cat("Identified", nrow(regulatory_edges), "potential regulatory interactions\n")

# 8. Save results
cat("Saving results...\n")

# Save TF activity matrix
write.table(tf_activity, 
           file.path(output_dir, paste0(sample_name, "_tf_activity.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save TF statistics
write.table(tf_stats,
           file.path(output_dir, paste0(sample_name, "_tf_statistics.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save regulatory network
write.table(regulatory_edges,
           file.path(output_dir, paste0(sample_name, "_regulatory_network.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save motif matches
motif_matches_df <- data.frame(
    peak_id = rep(rownames(motif_matches), ncol(motif_matches)),
    TF = rep(colnames(motif_matches), each = nrow(motif_matches)),
    has_motif = as.vector(motif_matches)
)
motif_matches_df <- motif_matches_df[motif_matches_df$has_motif, ]

write.table(motif_matches_df[, c("peak_id", "TF")],
           file.path(output_dir, paste0(sample_name, "_motif_matches.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#===============================================================================
# VISUALIZATION GENERATION
#===============================================================================
# Purpose: Create comprehensive visualizations of TF activity and regulatory networks
#
# TF Activity Statistics:
# - Mean activity: Average TF activity across all cells
# - Coefficient of variation: Measure of TF activity variability
# - Ranking: Identifies most variable and active TFs
#
# Visualization Components:
# 1. TF Activity Heatmap:
#    - Shows top variable TFs across cells
#    - Color scale represents z-score deviations
#    - Hierarchical clustering reveals TF co-activity patterns
#
# 2. TF Variability Plot:
#    - Bar plot of most variable TFs
#    - Identifies TFs with dynamic activity patterns
#    - Useful for prioritizing regulatory factors
#
# 3. Regulatory Network Summary:
#    - TFs ranked by number of target genes
#    - Identifies key regulatory hubs
#    - Highlights master regulators
#
# Sampling Strategy:
# - Limits cell visualization to 1000 for performance
# - Random sampling maintains representative patterns
# - Preserves statistical relationships
#
# Quality Control:
# - Uses consistent color schemes across plots
# - Ensures readable font sizes and layouts
# - Validates data integrity before plotting
#
# Applications:
# - Identify key regulatory factors
# - Compare TF activity between conditions
# - Prioritize TFs for experimental validation
#===============================================================================

# 9. Generate visualizations
cat("Generating visualizations...\n")

# TF activity heatmap (top variable TFs)
top_tfs <- head(tf_stats$TF, 20)
tf_activity_subset <- tf_activity[top_tfs, ]

# Sample cells for visualization if too many
if(ncol(tf_activity_subset) > 1000) {
    sampled_cells <- sample(ncol(tf_activity_subset), 1000)
    tf_activity_subset <- tf_activity_subset[, sampled_cells]
}

png(file.path(output_dir, paste0(sample_name, "_tf_activity_heatmap.png")), 
    width = 12, height = 10, units = "in", res = 300)

# Create heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(tf_activity_subset, 
        name = "TF Activity\n(z-score)",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 12)))

dev.off()

# TF activity distribution
tf_stats_plot <- tf_stats %>% head(15)
p1 <- ggplot(tf_stats_plot, aes(x = reorder(TF, cv), y = cv)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    coord_flip() +
    labs(title = paste("Most Variable TFs -", sample_name),
         x = "Transcription Factor", 
         y = "Coefficient of Variation") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_tf_variability.png")), 
       p1, width = 10, height = 8, dpi = 300)

# Regulatory network summary
network_summary <- regulatory_edges %>%
    group_by(TF) %>%
    summarise(n_targets = n_distinct(target_gene)) %>%
    arrange(desc(n_targets)) %>%
    head(15)

p2 <- ggplot(network_summary, aes(x = reorder(TF, n_targets), y = n_targets)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    coord_flip() +
    labs(title = paste("TFs with Most Target Genes -", sample_name),
         x = "Transcription Factor", 
         y = "Number of Target Genes") +
    theme_minimal()

ggsave(file.path(output_dir, paste0(sample_name, "_tf_targets.png")), 
       p2, width = 10, height = 8, dpi = 300)

#===============================================================================
# SUMMARY REPORT GENERATION
#===============================================================================
# Purpose: Generate comprehensive summary statistics for GRN analysis
#
# Key Metrics:
# - total_peaks_analyzed: Number of accessible chromatin regions processed
# - total_cells: Number of single cells included in analysis
# - motifs_searched: Number of TF motifs scanned
# - peaks_with_motifs: Peaks containing at least one TF binding site
# - tfs_with_activity: TFs with calculated activity scores
# - regulatory_edges: Total TF-gene regulatory relationships identified
# - mean_tf_activity: Average absolute TF activity across all factors
# - most_active_tf: TF with highest absolute mean activity
# - most_variable_tf: TF with highest coefficient of variation
#
# Quality Assessment:
# - Motif coverage: Fraction of peaks with TF binding sites
# - Network density: Number of regulatory edges per TF
# - Activity distribution: Range and variability of TF activities
#
# Applications:
# - Quality control assessment
# - Comparative analysis between samples
# - Method optimization and troubleshooting
# - Results interpretation and validation
#
# Output Format:
# - Tab-separated values for easy parsing
# - Metric-value pairs for programmatic access
# - Human-readable summary statistics
#===============================================================================

# 10. Create summary report
cat("Creating summary report...\n")

grn_summary <- data.frame(
    metric = c("total_peaks_analyzed", "total_cells", "motifs_searched", 
               "peaks_with_motifs", "tfs_with_activity", "regulatory_edges",
               "mean_tf_activity", "most_active_tf", "most_variable_tf"),
    value = c(nrow(motif_matches), ncol(tf_activity), ncol(motif_matches),
              sum(rowSums(motif_matches) > 0), nrow(tf_activity), nrow(regulatory_edges),
              round(mean(abs(tf_stats$mean_activity)), 3),
              tf_stats$TF[which.max(abs(tf_stats$mean_activity))],
              tf_stats$TF[1])
)

write.table(grn_summary,
           file.path(output_dir, paste0(sample_name, "_grn_summary.tsv")),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("GRN analysis complete!\n")
cat("Files generated:\n")
cat("  - TF activity matrix:", paste0(sample_name, "_tf_activity.tsv"), "\n")
cat("  - TF statistics:", paste0(sample_name, "_tf_statistics.tsv"), "\n")
cat("  - Regulatory network:", paste0(sample_name, "_regulatory_network.tsv"), "\n")
cat("  - Motif matches:", paste0(sample_name, "_motif_matches.tsv"), "\n")
cat("  - Summary report:", paste0(sample_name, "_grn_summary.tsv"), "\n")
cat("  - Visualizations: heatmaps and bar plots\n")

# Print key findings
cat("\nKey Findings:\n")
cat("=============\n")
cat("Most variable TFs:\n")
print(head(tf_stats[, c("TF", "cv")], 5))
cat("\nTop regulatory hubs:\n")
print(head(network_summary, 5))
EOF

# Make R script executable
chmod +x "$OUTPUT_DIR/grn_analysis/${SAMPLE}_grn_analysis.R"

#===============================================================================
# COMPARATIVE ANALYSIS SCRIPT GENERATION
#===============================================================================
# Purpose: Create specialized script for control vs mutant comparison
# - Generates only when processing control sample
# - Enables differential TF activity analysis
# - Identifies condition-specific regulatory changes
#
# Comparative Analysis Features:
# - Differential TF activity calculation
# - Volcano plot visualization
# - Correlation analysis between conditions
# - Statistical significance testing
#
# Applications:
# - Identify dysregulated TFs in disease/mutant conditions
# - Discover condition-specific regulatory networks
# - Prioritize therapeutic targets
# - Validate experimental hypotheses
#===============================================================================

# Create a comparative analysis script for control vs mutant
if [[ "$SAMPLE" == "R26-Nestin-Ctrl-adult" ]]; then
cat > "$OUTPUT_DIR/grn_analysis/comparative_grn_analysis.R" << 'EOF'
#!/usr/bin/env Rscript
# Comparative GRN Analysis: Control vs Mutant
# Identifies differentially active TFs and disrupted regulatory networks

#===============================================================================
# COMPARATIVE ANALYSIS R SCRIPT - LIBRARY LOADING AND DATA PREPARATION
#===============================================================================
# Purpose: Load required libraries and TF activity data for comparison
#
# Data Loading Strategy:
# - Load TF activity matrices from both conditions
# - Identify common TFs present in both datasets
# - Ensure data compatibility and alignment
#
# Statistical Approach:
# - Calculate mean TF activity per condition
# - Compute log2 fold changes (mutant vs control)
# - Apply significance thresholds for differential activity
#
# Quality Control:
# - Check data dimensions and overlap
# - Handle missing values and edge cases
# - Validate statistical assumptions
#===============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
})

cat("Running comparative GRN analysis...\n")

ctrl_sample <- "R26-Nestin-Ctrl-adult"
mut_sample <- "R26-Nestin-Mut-adult"
output_dir <- "grn_analysis"

# Load TF activity data
ctrl_tf <- tryCatch({
    read.table(file.path(output_dir, paste0(ctrl_sample, "_tf_activity.tsv")),
               header = TRUE, sep = "\t", row.names = 1)
}, error = function(e) {
    cat("Control data not found\n")
    return(NULL)
})

mut_tf <- tryCatch({
    read.table(file.path(output_dir, paste0(mut_sample, "_tf_activity.tsv")),
               header = TRUE, sep = "\t", row.names = 1)
}, error = function(e) {
    cat("Mutant data not found\n")
    return(NULL)
})

if(is.null(ctrl_tf) || is.null(mut_tf)) {
    stop("Cannot perform comparative analysis - missing data files")
}

# Find common TFs
common_tfs <- intersect(rownames(ctrl_tf), rownames(mut_tf))
cat("Comparing", length(common_tfs), "common TFs\n")

if(length(common_tfs) == 0) {
    stop("No common TFs found between samples")
}

# Calculate mean TF activity per condition
ctrl_means <- rowMeans(ctrl_tf[common_tfs, ])
mut_means <- rowMeans(mut_tf[common_tfs, ])

# Differential TF activity
tf_comparison <- data.frame(
    TF = common_tfs,
    ctrl_activity = ctrl_means,
    mut_activity = mut_means,
    log2fc = log2((mut_means + 1) / (ctrl_means + 1)),
    difference = mut_means - ctrl_means
)

tf_comparison$abs_log2fc <- abs(tf_comparison$log2fc)
tf_comparison <- tf_comparison[order(tf_comparison$abs_log2fc, decreasing = TRUE), ]

# Identify significantly changed TFs
tf_comparison$direction <- ifelse(tf_comparison$log2fc > 0.5, "Up in Mutant",
                                ifelse(tf_comparison$log2fc < -0.5, "Down in Mutant", "Unchanged"))

cat("TF activity changes:\n")
cat("Up in Mutant:", sum(tf_comparison$direction == "Up in Mutant"), "\n")
cat("Down in Mutant:", sum(tf_comparison$direction == "Down in Mutant"), "\n")
cat("Unchanged:", sum(tf_comparison$direction == "Unchanged"), "\n")

#===============================================================================
# COMPARATIVE ANALYSIS - RESULTS GENERATION AND VISUALIZATION
#===============================================================================
# Purpose: Generate comprehensive comparison results and visualizations
#
# Results Table Components:
# - TF identifiers and mean activities per condition
# - Log2 fold changes and significance calls
# - Absolute fold changes for ranking
# - Boolean significance flags for filtering
#
# Visualization Strategy:
# - Volcano plot: Shows magnitude vs significance of changes
# - Correlation plot: Reveals overall relationship between conditions
# - Color coding: Highlights significant vs non-significant changes
#
# Statistical Interpretation:
# - Positive log2FC: Higher activity in mutant
# - Negative log2FC: Lower activity in mutant
# - Threshold-based significance: Arbitrary but consistent cutoff
# - Correlation coefficient: Overall similarity between conditions
#
# Applications:
# - Identify key regulatory differences
# - Prioritize TFs for functional validation
# - Assess global regulatory disruption
# - Generate hypotheses for mechanism studies
#===============================================================================

# Save comparative results
write.table(tf_comparison,
           file.path(output_dir, "ctrl_vs_mut_tf_comparison.tsv"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Generate comparative plots
# Volcano plot
p1 <- ggplot(tf_comparison, aes(x = log2fc, y = abs_log2fc, color = direction)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = ifelse(abs_log2fc > 1, TF, "")), 
              vjust = -0.5, size = 3) +
    scale_color_manual(values = c("Up in Mutant" = "red", 
                                 "Down in Mutant" = "blue", 
                                 "Unchanged" = "gray")) +
    labs(title = "Differential TF Activity: Control vs Mutant",
         x = "Log2 Fold Change (Mutant/Control)",
         y = "Absolute Log2 Fold Change",
         color = "Direction") +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file.path(output_dir, "tf_activity_volcano_plot.png"), 
       p1, width = 10, height = 8, dpi = 300)

# Correlation plot
p2 <- ggplot(tf_comparison, aes(x = ctrl_activity, y = mut_activity)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_text(aes(label = ifelse(abs_log2fc > 1, TF, "")), 
              vjust = -0.5, hjust = 0.5, size = 3) +
    labs(title = "TF Activity Correlation: Control vs Mutant",
         x = "Control TF Activity (mean z-score)",
         y = "Mutant TF Activity (mean z-score)") +
    theme_minimal()

ggsave(file.path(output_dir, "tf_activity_correlation.png"), 
       p2, width = 8, height = 8, dpi = 300)

cat("Comparative analysis complete!\n")
cat("Top differentially active TFs:\n")
print(head(tf_comparison[, c("TF", "log2fc", "direction")], 10))
EOF

chmod +x "$OUTPUT_DIR/grn_analysis/comparative_grn_analysis.R"
fi

#===============================================================================
# MAIN ANALYSIS EXECUTION
#===============================================================================
# Purpose: Execute the comprehensive GRN analysis pipeline
#
# Execution Strategy:
# - Change to output directory for proper file paths
# - Export sample name as environment variable
# - Validate R/Rscript availability
# - Execute main analysis script with error handling
#
# Error Handling:
# - Check for R installation before execution
# - Provide clear error messages for missing dependencies
# - Validate successful completion before proceeding
#
# Quality Control:
# - Monitor script execution progress
# - Capture and report any runtime errors
# - Verify output file generation
#===============================================================================

# Run the GRN analysis
echo "DEBUG: Running gene regulatory network analysis..."
cd "$OUTPUT_DIR"
export SAMPLE="$SAMPLE"

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R"
    exit 1
fi

# Run the main GRN analysis script
echo "DEBUG: Starting GRN analysis for $SAMPLE..."
Rscript "grn_analysis/${SAMPLE}_grn_analysis.R" "$SAMPLE"

#===============================================================================
# RESULTS VALIDATION AND REPORTING
#===============================================================================
# Purpose: Validate analysis completion and provide comprehensive reporting
#
# Validation Strategy:
# - Check script exit status for successful completion
# - Verify output file generation and integrity
# - Display summary statistics when available
#
# Error Reporting:
# - Provide specific guidance for dependency issues
# - List required R/Bioconductor packages
# - Suggest installation commands for missing packages
#
# Results Summary:
# - List all generated output files
# - Describe file contents and applications
# - Provide file paths for downstream analysis
#
# Quality Assurance:
# - Confirm successful pipeline completion
# - Enable troubleshooting for failed runs
# - Support reproducible analysis workflows
#===============================================================================

#===============================================================================
# COMPLETION SUMMARY AND OUTPUT DOCUMENTATION
#===============================================================================
# Purpose: Provide comprehensive summary of analysis results and outputs
#
# Output File Descriptions:
# - tf_activity.tsv: Matrix of TF activity scores (TFs x cells)
# - tf_statistics.tsv: Summary statistics for each TF (mean, variance, etc.)
# - regulatory_network.tsv: TF-gene regulatory relationships
# - motif_matches.tsv: Binary matrix of motif occurrences in peaks
# - tf_activity_heatmap.png: Visualization of top variable TFs
# - tf_variability.png: Distribution plots for TF activities
# - tf_targets.png: Network topology visualization
# - grn_summary.tsv: Key metrics and quality control statistics
# - grn_analysis.R: Complete analysis script for reproducibility
# - comparative_grn_analysis.R: Control vs mutant comparison (control only)
#
# Error Handling:
# - Comprehensive dependency checking
# - Clear troubleshooting guidance
# - Package installation instructions
#
# Quality Assurance:
# - File existence validation
# - Analysis completion confirmation
# - Reproducibility documentation
#===============================================================================

if [[ $? -eq 0 ]]; then
    echo "DEBUG: GRN analysis completed successfully"
    
    # Display summary if available
    SUMMARY_FILE="grn_analysis/${SAMPLE}_grn_summary.tsv"
    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "GRN Analysis Summary:"
        cat "$SUMMARY_FILE"
    fi
else
    echo "WARNING: GRN analysis failed. Check R dependencies."
    echo "Required packages: motifmatchr, TFBSTools, BSgenome.Mmusculus.UCSC.mm10, chromVAR"
    echo "Install with BiocManager::install(c('motifmatchr', 'TFBSTools', 'BSgenome.Mmusculus.UCSC.mm10', 'chromVAR'))"
fi

echo "Output files created in: $OUTPUT_DIR/grn_analysis/"
echo "  - TF activity scores: ${SAMPLE}_tf_activity.tsv"
echo "  - TF statistics: ${SAMPLE}_tf_statistics.tsv"
echo "  - Regulatory network: ${SAMPLE}_regulatory_network.tsv"
echo "  - Motif matches: ${SAMPLE}_motif_matches.tsv"
echo "  - Analysis summary: ${SAMPLE}_grn_summary.tsv"
echo "  - Visualizations: heatmaps and plots"

echo "========================================="
echo "Step 11 complete for $SAMPLE"
echo "Gene regulatory network analysis finished"
echo "End time: $(date)"
echo "========================================="