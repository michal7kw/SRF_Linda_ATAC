#!/bin/bash

# Test zcat in conda environment
echo "=== Testing zcat in conda environment ==="

# Setup conda environment like the pipeline
if [[ -f "/opt/common/tools/ric.cosr/miniconda3/bin/activate" ]]; then
    source /opt/common/tools/ric.cosr/miniconda3/bin/activate
    conda activate peak_calling_new
    echo "✓ Conda environment 'peak_calling_new' activated"
else
    echo "✗ Conda activation script not found"
    exit 1
fi

# Test commands in conda environment
echo "Testing which zcat..."
which zcat

echo "Testing zcat version..."
zcat --version 2>&1 || echo "No version info available"

echo "Testing PATH:"
echo "$PATH"

# Test the actual file
FEATURES_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output/filtered_peak_bc_matrix/features.tsv.gz"
echo "Testing file: $FEATURES_FILE"

if [[ -f "$FEATURES_FILE" ]]; then
    echo "✓ File exists in conda environment"
    echo "Testing zcat decompression in conda environment..."
    if zcat "$FEATURES_FILE" | head -n 1 > /dev/null 2>&1; then
        echo "✓ zcat works in conda environment"
    else
        echo "✗ zcat fails in conda environment"
        echo "Exit code: $?"
        
        # Try alternative
        echo "Testing with gunzip -dc..."
        if gunzip -dc "$FEATURES_FILE" | head -n 1 > /dev/null 2>&1; then
            echo "✓ gunzip -dc works"
        else
            echo "✗ gunzip -dc also fails"
        fi
    fi
else
    echo "✗ File not found in conda environment"
fi