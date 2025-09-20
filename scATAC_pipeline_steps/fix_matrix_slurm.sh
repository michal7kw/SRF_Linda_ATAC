#!/bin/bash
#SBATCH --job-name=fix_matrix_format
#SBATCH --output=logs/fix_matrix_%j.out
#SBATCH --error=logs/fix_matrix_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# Fix the matrix.mtx format issue in SLURM environment
echo "[INFO $(date)] === Fixing Matrix Format Issue ==="

MATRIX_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix"
COUNTS_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/peak_barcode_counts.txt"
TEMP_DIR="/beegfs/scratch/tmp/fix_matrix_$$"

echo "[INFO $(date)] Working with matrix directory: $MATRIX_DIR"
echo "[INFO $(date)] Input counts file: $COUNTS_FILE"
echo "[INFO $(date)] Temporary directory: $TEMP_DIR"

mkdir -p "$TEMP_DIR"

# Extract existing files
echo "[INFO $(date)] Extracting features and barcodes..."
gunzip -dc "$MATRIX_DIR/features.tsv.gz" > "$TEMP_DIR/features.tsv"
gunzip -dc "$MATRIX_DIR/barcodes.tsv.gz" > "$TEMP_DIR/barcodes.tsv"

# Create index mappings
echo "[INFO $(date)] Creating index mappings..."
awk '{print $1 "\t" NR}' "$TEMP_DIR/features.tsv" > "$TEMP_DIR/feature_index.txt"
awk '{print $1 "\t" NR}' "$TEMP_DIR/barcodes.tsv" > "$TEMP_DIR/barcode_index.txt"

# Get dimensions
num_features=$(wc -l < "$TEMP_DIR/features.tsv")
num_barcodes=$(wc -l < "$TEMP_DIR/barcodes.tsv")
num_entries=$(wc -l < "$COUNTS_FILE")

echo "[INFO $(date)] Matrix dimensions: $num_features features × $num_barcodes barcodes"
echo "[INFO $(date)] Number of entries: $num_entries"

# Create corrected matrix
echo "[INFO $(date)] Creating corrected matrix.mtx..."
{
    echo "%%MatrixMarket matrix coordinate integer general"
    echo "%"
    echo "$num_features $num_barcodes $num_entries"
    
    # Process counts file with correct joins and field mapping
    # After joins: barcode_pos peak_coord count feature_idx barcode_idx
    join -1 1 -2 1 <(sort -k1,1 "$COUNTS_FILE") <(sort -k1,1 "$TEMP_DIR/feature_index.txt") | \
    join -1 2 -2 1 <(sort -k2,2) <(sort -k1,1 "$TEMP_DIR/barcode_index.txt") | \
    awk '{print $4, $5, $3}' | \
    sort -k1,1n -k2,2n
    
} > "$TEMP_DIR/matrix_fixed.mtx"

# Verify the format
echo "[INFO $(date)] Verifying matrix format..."
echo "First 5 entries:"
tail -n +4 "$TEMP_DIR/matrix_fixed.mtx" | head -n 5

# Check format of sample
invalid_entries=$(tail -n +4 "$TEMP_DIR/matrix_fixed.mtx" | head -n 1000 | awk 'NF != 3 || $1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= 0 {count++} END {print count+0}')

if [[ $invalid_entries -gt 0 ]]; then
    echo "[ERROR $(date)] Still have $invalid_entries invalid entries"
    exit 1
fi

echo "[INFO $(date)] Matrix format is now correct!"

# Replace the compressed matrix file
echo "[INFO $(date)] Compressing and replacing the matrix file..."
gzip "$TEMP_DIR/matrix_fixed.mtx"
mv "$TEMP_DIR/matrix_fixed.mtx.gz" "$MATRIX_DIR/matrix.mtx.gz"

echo "[INFO $(date)] Matrix file has been fixed!"

# Cleanup
echo "[INFO $(date)] Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "[INFO $(date)] Fix completed successfully!"