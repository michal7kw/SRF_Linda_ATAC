#!/bin/bash

# Fix the matrix.mtx format issue
echo "=== Fixing Matrix Format Issue ==="

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps

MATRIX_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix"
COUNTS_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/peak_barcode_counts.txt"

echo "Working with matrix directory: $MATRIX_DIR"
echo "Input counts file: $COUNTS_FILE"

# Extract existing files
echo "Extracting features and barcodes..."
gunzip -dc "$MATRIX_DIR/features.tsv.gz" > /tmp/features.tsv
gunzip -dc "$MATRIX_DIR/barcodes.tsv.gz" > /tmp/barcodes.tsv

# Create index mappings
echo "Creating index mappings..."
awk '{print $1 "\t" NR}' /tmp/features.tsv > /tmp/feature_index.txt
awk '{print $1 "\t" NR}' /tmp/barcodes.tsv > /tmp/barcode_index.txt

# Get dimensions
num_features=$(wc -l < /tmp/features.tsv)
num_barcodes=$(wc -l < /tmp/barcodes.tsv)
num_entries=$(wc -l < "$COUNTS_FILE")

echo "Matrix dimensions: $num_features features × $num_barcodes barcodes"
echo "Number of entries: $num_entries"

# Create corrected matrix
echo "Creating corrected matrix.mtx..."
{
    echo "%%MatrixMarket matrix coordinate integer general"
    echo "%"
    echo "$num_features $num_barcodes $num_entries"
    
    # Process counts file with correct joins
    join -1 1 -2 1 <(sort -k1,1 "$COUNTS_FILE") <(sort -k1,1 /tmp/feature_index.txt) | \
    join -1 2 -2 1 <(sort -k2,2) <(sort -k1,1 /tmp/barcode_index.txt) | \
    awk '{print $4, $3, $2}' | \
    sort -k1,1n -k2,2n
    
} > /tmp/matrix_fixed.mtx

# Verify the format
echo "Verifying matrix format..."
echo "First 5 entries:"
tail -n +4 /tmp/matrix_fixed.mtx | head -n 5

# Check format
invalid_entries=$(tail -n +4 /tmp/matrix_fixed.mtx | head -n 1000 | awk 'NF != 3 || $1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= 0 {count++} END {print count+0}')

if [[ $invalid_entries -gt 0 ]]; then
    echo "ERROR: Still have $invalid_entries invalid entries"
    exit 1
fi

echo "Matrix format is now correct!"

# Replace the compressed matrix file
echo "Replacing the matrix file..."
gzip /tmp/matrix_fixed.mtx
mv /tmp/matrix_fixed.mtx.gz "$MATRIX_DIR/matrix.mtx.gz"

echo "Matrix file has been fixed!"

# Cleanup
rm -f /tmp/features.tsv /tmp/barcodes.tsv /tmp/feature_index.txt /tmp/barcode_index.txt

echo "Fix completed successfully!"