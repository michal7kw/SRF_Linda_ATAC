#!/bin/bash

# Debug script to test file validation
echo "=== Debug Validation Script ==="
echo "Working directory: $(pwd)"
echo "Current user: $(whoami)"

# Test file existence
FEATURES_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix/features.tsv.gz"
echo "Testing file: $FEATURES_FILE"

if [[ -f "$FEATURES_FILE" ]]; then
    echo "✓ File exists"
    echo "File size: $(stat -c%s "$FEATURES_FILE") bytes"
    echo "File permissions: $(stat -c%A "$FEATURES_FILE")"
    echo "File owner: $(stat -c%U:%G "$FEATURES_FILE")"
    
    # Test compression
    echo "Testing zcat decompression..."
    if zcat "$FEATURES_FILE" | head -n 1 > /dev/null 2>&1; then
        echo "✓ File can be decompressed"
        echo "First line:"
        zcat "$FEATURES_FILE" | head -n 1
    else
        echo "✗ File cannot be decompressed"
        echo "Exit code: $?"
        
        # Try alternative decompression methods
        echo "Testing gunzip -t (test integrity)..."
        gunzip -t "$FEATURES_FILE" 2>&1 || echo "Integrity test failed"
        
        echo "Testing file command..."
        file "$FEATURES_FILE"
        
        echo "Testing hexdump of first 32 bytes..."
        hexdump -C "$FEATURES_FILE" | head -n 2
    fi
else
    echo "✗ File does not exist"
    echo "Directory listing:"
    ls -la "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix/"
fi