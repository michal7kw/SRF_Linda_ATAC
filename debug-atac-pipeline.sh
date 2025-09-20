#!/bin/bash
#SBATCH --job-name=atac_debug
#SBATCH --output=logs/atac_debug_%A.out
#SBATCH --error=logs/atac_debug_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

# Robust ATAC processing pipeline with debugging
# Remove 'set -e' to continue on errors and see where problems occur

# Configuration
SAMPLE="R26-Nestin-Mut-adult"  # Modify as needed
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/quality_filtered_output_debug"
REF_GENOME="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/mm10/mm10.fa"
CHROMAP_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/mm10/chromap_index/genome.index"
THREADS=16

# Create directories first
echo "Creating directories..."
mkdir -p "$OUTPUT_DIR/logs" "$OUTPUT_DIR/qc" || {
    echo "ERROR: Failed to create directories"
    exit 1
}

echo "========================================="
echo "Processing sample: $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Step 1: Assess data quality with error checking
echo "Step 1: Assessing data quality..."

# First, check if the file exists and is readable
if [ ! -f "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" ]; then
    echo "ERROR: Input file not found: $DATA_DIR/${SAMPLE}_R2_001.fastq.gz"
    exit 1
fi

# Check file size
FILE_SIZE=$(stat -c%s "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz")
echo "R2 file size: $FILE_SIZE bytes"

# Try extracting a small number first
echo "Testing with 100 barcodes first..."
if zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" 2>/dev/null | \
   awk 'NR%4==2 {print substr($0,1,16)}' | \
   head -100 > "$OUTPUT_DIR/qc/test_100.txt"; then
    echo "Successfully extracted 100 barcodes"
    head -5 "$OUTPUT_DIR/qc/test_100.txt"
else
    echo "ERROR: Failed to extract test barcodes"
    echo "Trying alternative method..."
    gunzip -c "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" 2>/dev/null | \
        sed -n '2~4p' | cut -c1-16 | head -100 > "$OUTPUT_DIR/qc/test_100.txt"
fi

# Now try the full 100k extraction with timeout
echo "Extracting 100k barcodes (with 5 minute timeout)..."
timeout 300 bash -c "zcat '$DATA_DIR/${SAMPLE}_R2_001.fastq.gz' | \
    awk 'NR%4==2 {print substr(\$0,1,16)}' | \
    head -100000 > '$OUTPUT_DIR/qc/barcodes_100k.txt'"

if [ $? -eq 124 ]; then
    echo "WARNING: Barcode extraction timed out. Using smaller sample..."
    timeout 60 bash -c "zcat '$DATA_DIR/${SAMPLE}_R2_001.fastq.gz' | \
        awk 'NR%4==2 {print substr(\$0,1,16)}' | \
        head -10000 > '$OUTPUT_DIR/qc/barcodes_100k.txt'"
fi

# Check if file was created
if [ ! -f "$OUTPUT_DIR/qc/barcodes_100k.txt" ] || [ ! -s "$OUTPUT_DIR/qc/barcodes_100k.txt" ]; then
    echo "ERROR: Failed to create barcode file. Using fallback method..."
    # Fallback: use sed instead of awk
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" 2>/dev/null | \
        sed -n '2~4p' | cut -c1-16 | \
        head -10000 > "$OUTPUT_DIR/qc/barcodes_100k.txt" || {
            echo "FATAL: Cannot extract barcodes from R2 file"
            exit 1
        }
fi

# Count N bases
TOTAL_BC=$(wc -l < "$OUTPUT_DIR/qc/barcodes_100k.txt")
echo "Total barcodes extracted: $TOTAL_BC"

if [ $TOTAL_BC -eq 0 ]; then
    echo "ERROR: No barcodes extracted. Check input file."
    exit 1
fi

N_BC=$(grep -c 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" || echo "0")
CLEAN_BC=$(grep -vc 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" || echo "0")

# Handle division by zero
if [ $TOTAL_BC -gt 0 ]; then
    N_PERCENT=$(echo "scale=2; $N_BC * 100 / $TOTAL_BC" | bc)
else
    N_PERCENT=0
fi

echo "Barcode quality summary:"
echo "  Total barcodes: $TOTAL_BC"
echo "  Barcodes with N's: $N_BC ($N_PERCENT%)"
echo "  Clean barcodes: $CLEAN_BC"

# Get whitelist files
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
mkdir -p "$WHITELIST_DIR"

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1.txt" ]]; then
    echo "Downloading barcode whitelist..."
    wget -O "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz || {
            echo "ERROR: Failed to download whitelist"
            exit 1
        }
    gunzip "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
fi

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

# Check orientation using only clean barcodes
echo "Checking whitelist orientation..."
if [ $CLEAN_BC -gt 0 ]; then
    grep -v 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" | head -10000 > "$OUTPUT_DIR/qc/clean_barcodes_10k.txt" || {
        echo "WARNING: Could not extract clean barcodes"
        head -10000 "$OUTPUT_DIR/qc/barcodes_100k.txt" > "$OUTPUT_DIR/qc/clean_barcodes_10k.txt"
    }
else
    echo "WARNING: No clean barcodes found. Using all barcodes..."
    head -10000 "$OUTPUT_DIR/qc/barcodes_100k.txt" > "$OUTPUT_DIR/qc/clean_barcodes_10k.txt"
fi

MATCHES_ORIG=$(grep -Ff "$WHITELIST_DIR/atac_737K-arc-v1.txt" "$OUTPUT_DIR/qc/clean_barcodes_10k.txt" | wc -l || echo "0")
MATCHES_RC=$(grep -Ff "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" "$OUTPUT_DIR/qc/clean_barcodes_10k.txt" | wc -l || echo "0")

echo "Barcode whitelist matching:"
echo "  Original whitelist: $MATCHES_ORIG / $(wc -l < $OUTPUT_DIR/qc/clean_barcodes_10k.txt)"
echo "  Reverse complement: $MATCHES_RC / $(wc -l < $OUTPUT_DIR/qc/clean_barcodes_10k.txt)"

if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    echo "Using original whitelist"
fi

# Quick chromap test without quality filtering
echo ""
echo "Step 2: Running quick chromap test..."

# Extract 16bp barcodes
echo "Extracting 16bp barcodes..."
zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" 2>/dev/null | \
    awk 'NR%4==1 {print $0} 
         NR%4==2 {print substr($0,1,16)} 
         NR%4==3 {print "+"} 
         NR%4==0 {print substr($0,1,16)}' | \
    head -4000000 | \
    gzip > "$OUTPUT_DIR/${SAMPLE}_R2_16bp_test.fastq.gz" || {
        echo "ERROR: Failed to extract barcodes"
        exit 1
    }

# Test with chromap on subset
echo "Running chromap on subset..."
if command -v chromap &> /dev/null; then
    # Create subset of R1 and R3 for testing
    zcat "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" | head -4000000 | gzip > "$OUTPUT_DIR/${SAMPLE}_R1_test.fastq.gz"
    zcat "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" | head -4000000 | gzip > "$OUTPUT_DIR/${SAMPLE}_R3_test.fastq.gz"
    
    chromap --preset atac \
        -x "$CHROMAP_INDEX" \
        -r "$REF_GENOME" \
        -1 "$OUTPUT_DIR/${SAMPLE}_R1_test.fastq.gz" \
        -2 "$OUTPUT_DIR/${SAMPLE}_R3_test.fastq.gz" \
        -b "$OUTPUT_DIR/${SAMPLE}_R2_16bp_test.fastq.gz" \
        --barcode-whitelist "$WHITELIST" \
        --barcode-mismatch 2 \
        --read-mismatch 5 \
        --threads $THREADS \
        -o "$OUTPUT_DIR/${SAMPLE}_fragments_test.tsv" \
        --low-mem 2>&1 | tee "$OUTPUT_DIR/chromap_test.log"
    
    if [ -f "$OUTPUT_DIR/${SAMPLE}_fragments_test.tsv" ]; then
        FRAGS=$(wc -l < "$OUTPUT_DIR/${SAMPLE}_fragments_test.tsv")
        echo "Test fragments generated: $FRAGS"
        
        # Check valid barcode rate
        VALID_BC_TEST=$(cut -f4 "$OUTPUT_DIR/${SAMPLE}_fragments_test.tsv" | \
                        grep -Ff "$WHITELIST" | wc -l || echo "0")
        TOTAL_BC_TEST=$(cut -f4 "$OUTPUT_DIR/${SAMPLE}_fragments_test.tsv" | wc -l || echo "1")
        
        if [ $TOTAL_BC_TEST -gt 0 ]; then
            VALID_RATE=$(echo "scale=2; $VALID_BC_TEST * 100 / $TOTAL_BC_TEST" | bc)
            echo "Valid barcode rate in test: $VALID_RATE%"
        fi
    else
        echo "WARNING: Chromap test failed to generate fragments"
    fi
else
    echo "WARNING: chromap not found. Install from: https://github.com/haowenz/chromap"
fi

echo ""
echo "========================================="
echo "Diagnostic Summary:"
echo "  Input file exists: Yes"
echo "  Barcodes extracted: $TOTAL_BC"
echo "  N base rate: $N_PERCENT%"
echo "  Best whitelist: $(basename $WHITELIST)"
echo "  Match rate: $(echo "scale=2; $MATCHES_ORIG * 100 / 10000" | bc)%"
echo "========================================="

echo "Debug information saved to: $OUTPUT_DIR"
echo "End time: $(date)"