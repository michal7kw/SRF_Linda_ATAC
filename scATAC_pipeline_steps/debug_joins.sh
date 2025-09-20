#!/bin/bash

# Debug the join operations with small sample
echo "=== Debugging Join Operations ==="

MATRIX_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix"
COUNTS_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/peak_barcode_counts.txt"

# Create small test files
echo "Creating test files..."
head -n 5 "$COUNTS_FILE" > /tmp/test_counts.txt
gunzip -dc "$MATRIX_DIR/features.tsv.gz" | head -n 5 > /tmp/test_features.tsv
gunzip -dc "$MATRIX_DIR/barcodes.tsv.gz" | head -n 5 > /tmp/test_barcodes.tsv

echo "Test counts file:"
cat /tmp/test_counts.txt
echo ""

echo "Test features file:"
cat /tmp/test_features.tsv
echo ""

echo "Test barcodes file:"
cat /tmp/test_barcodes.tsv
echo ""

# Create index files
echo "Creating index mappings..."
awk '{print $1 "\t" NR}' /tmp/test_features.tsv > /tmp/test_feature_index.txt
awk '{print $1 "\t" NR}' /tmp/test_barcodes.tsv > /tmp/test_barcode_index.txt

echo "Feature index mapping:"
cat /tmp/test_feature_index.txt
echo ""

echo "Barcode index mapping:"
cat /tmp/test_barcode_index.txt
echo ""

# Test first join
echo "First join (counts + feature indices):"
join -1 1 -2 1 <(sort -k1,1 /tmp/test_counts.txt) <(sort -k1,1 /tmp/test_feature_index.txt)
echo ""

# Test second join
echo "Second join (result + barcode indices):"
join -1 1 -2 1 <(sort -k1,1 /tmp/test_counts.txt) <(sort -k1,1 /tmp/test_feature_index.txt) | \
join -1 2 -2 1 <(sort -k2,2) <(sort -k1,1 /tmp/test_barcode_index.txt)
echo ""

# Test final awk
echo "Final result with corrected awk:"
join -1 1 -2 1 <(sort -k1,1 /tmp/test_counts.txt) <(sort -k1,1 /tmp/test_feature_index.txt) | \
join -1 2 -2 1 <(sort -k2,2) <(sort -k1,1 /tmp/test_barcode_index.txt) | \
awk '{print "Fields:", NF, "| $1=" $1, "$2=" $2, "$3=" $3, "$4=" $4, "$5=" $5}'

# Cleanup
rm -f /tmp/test_*.txt /tmp/test_*.tsv