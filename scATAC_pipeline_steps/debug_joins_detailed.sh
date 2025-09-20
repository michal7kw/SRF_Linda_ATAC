#!/bin/bash

# More detailed debug of join operations
echo "=== Detailed Join Debug ==="

COUNTS_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/peak_barcode_counts.txt"
MATRIX_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix"

# Get actual barcodes from real files
echo "Getting barcodes from counts file vs barcodes.tsv:"
head -n 3 "$COUNTS_FILE" | cut -f2
echo "---"
gunzip -dc "$MATRIX_DIR/barcodes.tsv.gz" | head -n 3
echo ""

# Create small test with ACTUAL matching data
echo "Creating test files with matching barcodes..."
head -n 5 "$COUNTS_FILE" > /tmp/test_counts.txt
cut -f2 /tmp/test_counts.txt | sort -u > /tmp/test_barcodes.tsv
gunzip -dc "$MATRIX_DIR/features.tsv.gz" | head -n 5 > /tmp/test_features.tsv

echo "Test files created:"
echo "Counts:"
cat /tmp/test_counts.txt
echo ""
echo "Barcodes:"
cat /tmp/test_barcodes.tsv
echo ""

# Create indices
awk '{print $1 "\t" NR}' /tmp/test_features.tsv > /tmp/test_feature_index.txt
awk '{print $1 "\t" NR}' /tmp/test_barcodes.tsv > /tmp/test_barcode_index.txt

echo "Feature indices:"
cat /tmp/test_feature_index.txt
echo ""
echo "Barcode indices:"
cat /tmp/test_barcode_index.txt
echo ""

# Test the full pipeline
echo "=== Testing Full Join Pipeline ==="

echo "Step 1 - First join result:"
join -1 1 -2 1 <(sort -k1,1 /tmp/test_counts.txt) <(sort -k1,1 /tmp/test_feature_index.txt) > /tmp/step1.txt
cat /tmp/step1.txt
echo ""

echo "Step 2 - Second join input (sorted by field 2):"
sort -k2,2 /tmp/step1.txt
echo ""

echo "Step 2 - Second join result:"
join -1 2 -2 1 <(sort -k2,2 /tmp/step1.txt) <(sort -k1,1 /tmp/test_barcode_index.txt) > /tmp/step2.txt
cat /tmp/step2.txt
echo ""

echo "Step 3 - Final awk with field analysis:"
cat /tmp/step2.txt | awk '{print "NF=" NF " | $1=" $1 " $2=" $2 " $3=" $3 " $4=" $4 " $5=" $5 " | MTX format: " $4 " " $5 " " $3}'

# Cleanup
rm -f /tmp/test_*.txt /tmp/test_*.tsv /tmp/step*.txt