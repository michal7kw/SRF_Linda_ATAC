
# %%
OUTPUT_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

SAMPLES <- c("R26-Nestin-Ctrl-adult", "R26-Nestin-Mut-adult")
SAMPLE <- SAMPLES[1]

# %%
suppressPackageStartupMessages({
    library(Matrix)
    library(irlba)
    library(uwot)
    library(Rtsne)
    library(cluster)
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
})

sample_name <- SAMPLE

cat("Processing sample:", sample_name, "\n")

# Set paths
matrix_dir <- file.path(OUTPUT_DIR, "filtered_peak_bc_matrix")
output_dir <- file.path(OUTPUT_DIR, "dimensionality_reduction")

# Read 10X format data
cat("Reading peak-barcode matrix...\n")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
barcodes_path <- file.path(matrix_dir, "barcodes.tsv.gz")

# Check if files exist
if(!file.exists(matrix_path)) {
    stop("Matrix file not found: ", matrix_path)
}

# Read matrix
peak_matrix <- readMM(matrix_path)
features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Set names
rownames(peak_matrix) <- features[,1]  # Peak coordinates
colnames(peak_matrix) <- barcodes[,1]  # Cell barcodes

cat("Matrix dimensions:", nrow(peak_matrix), "peaks x", ncol(peak_matrix), "cells\n")

# Basic filtering
cat("Performing basic filtering...\n")

# Filter cells with too few or too many peaks
cell_peak_counts <- Matrix::colSums(peak_matrix > 0)
peak_cell_counts <- Matrix::rowSums(peak_matrix > 0)

# Debug: Show distribution of peak counts
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

# Use very conservative filtering based on actual data
min_peaks_per_cell <- max(1, quantile(cell_peak_counts[cell_peak_counts > 0], 0.01))
max_peaks_per_cell <- max(quantile(cell_peak_counts, 0.99), min_peaks_per_cell + 10)
min_cells_per_peak <- 1  # Very permissive

cat("Adjusted filtering thresholds:\n")
cat("  Min peaks per cell:", min_peaks_per_cell, "\n")
cat("  Max peaks per cell:", max_peaks_per_cell, "\n")
cat("  Min cells per peak:", min_cells_per_peak, "\n")

# Apply filters
good_cells <- which(cell_peak_counts >= min_peaks_per_cell & 
                   cell_peak_counts <= max_peaks_per_cell)
good_peaks <- which(peak_cell_counts >= min_cells_per_peak)

filtered_matrix <- peak_matrix[good_peaks, good_cells]
cat("After filtering:", nrow(filtered_matrix), "peaks x", ncol(filtered_matrix), "cells\n")

# TF-IDF normalization (LSI preprocessing)
cat("Performing TF-IDF normalization...\n")

# Convert to CsparseMatrix format if needed for TF-IDF calculation
if (class(filtered_matrix)[1] != "dgCMatrix") {
    filtered_matrix <- as(filtered_matrix, "CsparseMatrix")
}

# Calculate term frequency (TF)
tf_matrix <- filtered_matrix
tf_matrix@x <- tf_matrix@x / rep.int(Matrix::colSums(tf_matrix), diff(tf_matrix@p))

# Calculate inverse document frequency (IDF)
idf_values <- log(ncol(tf_matrix) / Matrix::rowSums(tf_matrix > 0))
idf_values[is.infinite(idf_values)] <- 0

# Apply IDF weighting
tfidf_matrix <- tf_matrix
for(i in 1:nrow(tfidf_matrix)) {
    tfidf_matrix[i, ] <- tfidf_matrix[i, ] * idf_values[i]
}

# LSI (Latent Semantic Indexing)
cat("Performing LSI...\n")
min_dim <- min(nrow(tfidf_matrix), ncol(tfidf_matrix))
max_safe_components <- min_dim - 3  # Even more conservative

# Start with a reasonable number and ensure it's safe
n_components <- min(50, max_safe_components)
n_components <- max(2, n_components)

# Final safety check
if (n_components >= min_dim) {
    n_components <- max(2, min_dim - 3)
}

cat("Matrix dimensions:", nrow(tfidf_matrix), "x", ncol(tfidf_matrix), "\n")
cat("Min dimension:", min_dim, ", Max safe components:", max_safe_components, "\n")
cat("Using", n_components, "LSI components\n")

# Additional check before calling irlba
if (n_components >= min_dim || n_components < 1) {
    stop("Invalid n_components: ", n_components, " for matrix of size ", 
         nrow(tfidf_matrix), "x", ncol(tfidf_matrix))
}

# Perform SVD using irlba for efficiency
set.seed(42)
lsi_result <- irlba(tfidf_matrix, nv = n_components, nu = n_components)

# LSI components (excluding first component which often correlates with sequencing depth)
lsi_embeddings <- lsi_result$v[, 2:n_components]
rownames(lsi_embeddings) <- colnames(filtered_matrix)
colnames(lsi_embeddings) <- paste0("LSI_", 1:(n_components-1))

cat("LSI completed. Using", ncol(lsi_embeddings), "components\n")

# PCA on LSI components
cat("Performing PCA on LSI components...\n")
pca_result <- prcomp(lsi_embeddings, center = TRUE, scale. = TRUE)
pca_embeddings <- pca_result$x[, 1:min(20, ncol(pca_result$x))]

# UMAP
cat("Computing UMAP embedding...\n")
set.seed(42)
umap_result <- umap(lsi_embeddings[, 1:min(30, ncol(lsi_embeddings))], 
                   n_neighbors = 30, 
                   min_dist = 0.3,
                   metric = "cosine",
                   verbose = TRUE)

umap_df <- data.frame(
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2],
    Barcode = rownames(lsi_embeddings)
)

# t-SNE (optional, commented out for speed)
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

# Clustering using k-means on LSI components
cat("Performing clustering...\n")
n_clusters <- 8  # Adjust based on expected cell types

set.seed(42)
cluster_result <- kmeans(lsi_embeddings[, 1:min(20, ncol(lsi_embeddings))], 
                        centers = n_clusters, 
                        nstart = 20)

umap_df$Cluster <- factor(cluster_result$cluster)

# Calculate silhouette score
sil_score <- silhouette(cluster_result$cluster, 
                       dist(lsi_embeddings[, 1:min(10, ncol(lsi_embeddings))]))
avg_sil_width <- mean(sil_score[, "sil_width"])
cat("Average silhouette width:", round(avg_sil_width, 3), "\n")

# Generate plots
cat("Generating plots...\n")

# UMAP plot colored by cluster
p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    labs(title = paste("UMAP Clustering -", sample_name),
         subtitle = paste("n =", nrow(umap_df), "cells,", n_clusters, "clusters")) +
    theme_minimal() +
    theme(legend.position = "right")

ggsave(file.path(output_dir, paste0(sample_name, "_umap_clusters.png")), 
       p1, width = 10, height = 8, dpi = 300)

# PCA variance explained plot
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2))[1:20]
pca_var_df <- data.frame(
    PC = 1:length(var_explained),
    Variance_Explained = var_explained,
    Cumulative = cumsum(var_explained)
)

p2 <- ggplot(pca_var_df, aes(x = PC, y = Variance_Explained)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_line(aes(y = Cumulative), color = "red", size = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(title = paste("PCA Variance Explained -", sample_name),
         x = "Principal Component", 
         y = "Proportion of Variance Explained") +
    theme_minimal() +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance"))

ggsave(file.path(output_dir, paste0(sample_name, "_pca_variance.png")), 
       p2, width = 10, height = 6, dpi = 300)

# Save results
cat("Saving results...\n")

# Save embeddings
write.table(lsi_embeddings, 
           file.path(output_dir, paste0(sample_name, "_lsi_embeddings.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(pca_embeddings, 
           file.path(output_dir, paste0(sample_name, "_pca_embeddings.tsv")),
           sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

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


# Make R script executable
chmod +x "$OUTPUT_DIR/dimensionality_reduction/${SAMPLE}_dim_reduction.R"