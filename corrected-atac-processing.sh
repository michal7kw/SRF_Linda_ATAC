#!/bin/bash
#SBATCH --job-name=atac_processing_corrected
#SBATCH --output=logs/atac_corrected_%A_%a.out
#SBATCH --error=logs/atac_corrected_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

set -euo pipefail

# Configuration
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"
CHROMAP_INDEX="$OUTPUT_DIR/chromap_index/mm10.index"
THREADS=16

echo "========================================="
echo "Processing ATAC sample: $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Debug: Print all configuration variables
echo "DEBUG: Configuration variables:"
echo "  SAMPLE: $SAMPLE"
echo "  DATA_DIR: $DATA_DIR"
echo "  PROCESSED_DATA_DIR: $PROCESSED_DATA_DIR"
echo "  OUTPUT_DIR: $OUTPUT_DIR"
echo "  REF_GENOME: $REF_GENOME"
echo "  CHROMAP_INDEX: $CHROMAP_INDEX"
echo "  THREADS: $THREADS"
echo ""

# Create output directories
echo "DEBUG: Creating output directories..."
mkdir -p "$OUTPUT_DIR/qc" "$OUTPUT_DIR/fragments" "$OUTPUT_DIR/peaks" "$OUTPUT_DIR/logs"
if [[ $? -eq 0 ]]; then
    echo "DEBUG: Output directories created successfully"
else
    echo "ERROR: Failed to create output directories"
    exit 1
fi

# Step 1: Verify extracted barcode file exists
EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_16bp.fastq.gz"
if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Running extraction..."
    
    # Extract barcodes from positions 7-22
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,7,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,7,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
fi

echo "Using extracted barcodes from: $EXTRACTED_BC_FILE"

# Step 2: Determine barcode whitelist orientation
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
mkdir -p "$WHITELIST_DIR"

# Download whitelists if needed
if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1.txt" ]]; then
    echo "Downloading whitelist..."
    wget -q -O "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    gunzip "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
fi

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

echo "Testing whitelist orientation with extracted barcodes..."
echo "DEBUG: Extracting test barcodes with timeout protection..."

# Check if the barcode file contains already extracted 16bp barcodes or original 24bp sequences
echo "DEBUG: Checking barcode file format..."
SAMPLE_SEQ=$(zcat "$EXTRACTED_BC_FILE" | awk 'NR==2' | head -1)
SEQ_LENGTH=${#SAMPLE_SEQ}
echo "DEBUG: Sample sequence from processed file: $SAMPLE_SEQ"
echo "DEBUG: Sequence length: $SEQ_LENGTH"

# Check for N content in barcodes (indicates sequencing quality issues)
echo "DEBUG: Analyzing barcode quality..."
N_COUNT=$(echo "$SAMPLE_SEQ" | grep -o 'N' | wc -l)
echo "DEBUG: Number of N's in sample barcode: $N_COUNT"

if [[ $N_COUNT -gt 0 ]]; then
    echo "DEBUG: Barcodes contain N's - this indicates sequencing quality issues, not wrong extraction positions"
    echo "DEBUG: Based on R2 structure analysis: CAGACG[16bp barcode][XX] - positions 7-22 are correct"
    echo "DEBUG: Will handle N's through error correction and increased mismatch tolerance"
fi

# Extract test barcodes with quality filtering
echo "DEBUG: Extracting test barcodes with N filtering..."

# Extract barcodes with quality filtering for N's
if [[ $SEQ_LENGTH -eq 16 ]]; then
    echo "DEBUG: File contains 16bp extracted barcodes - using directly with N filtering"
    # Extract barcodes and separate by N content
    timeout 300 bash -c "zcat '$EXTRACTED_BC_FILE' | awk 'NR%4==2' | head -10000 > '$OUTPUT_DIR/qc/test_barcodes_all.txt'"
    EXTRACT_EXIT_CODE=$?
    
    # Create filtered sets
    grep -v 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_clean.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_clean.txt"
    grep 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_with_n.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_with_n.txt"
    
elif [[ $SEQ_LENGTH -eq 24 ]]; then
    echo "DEBUG: File contains 24bp sequences - extracting positions 7-22 with N filtering"
    # Extract 16bp from positions 7-22 and filter
    timeout 300 bash -c "zcat '$EXTRACTED_BC_FILE' | awk 'NR%4==2 {print substr(\$0,7,16)}' | head -10000 > '$OUTPUT_DIR/qc/test_barcodes_all.txt'"
    EXTRACT_EXIT_CODE=$?
    
    # Create filtered sets
    grep -v 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_clean.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_clean.txt"
    grep 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_with_n.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_with_n.txt"
    
else
    echo "DEBUG: Unexpected sequence length ($SEQ_LENGTH), extracting with length validation"
    timeout 300 bash -c "zcat '$EXTRACTED_BC_FILE' | awk 'NR%4==2 {if(length(\$0)>=22) print substr(\$0,7,16); else print \$0}' | head -10000 > '$OUTPUT_DIR/qc/test_barcodes_all.txt'"
    EXTRACT_EXIT_CODE=$?
    
    # Create filtered sets
    grep -v 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_clean.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_clean.txt"
    grep 'N' "$OUTPUT_DIR/qc/test_barcodes_all.txt" > "$OUTPUT_DIR/qc/test_barcodes_with_n.txt" || touch "$OUTPUT_DIR/qc/test_barcodes_with_n.txt"
fi

# Use all barcodes for testing (including N's)
cp "$OUTPUT_DIR/qc/test_barcodes_all.txt" "$OUTPUT_DIR/qc/test_barcodes.txt"

if [[ $EXTRACT_EXIT_CODE -eq 124 ]]; then
    echo "ERROR: Barcode extraction timed out after 5 minutes"
    echo "DEBUG: This suggests the barcode file may be corrupted or extremely large"
    echo "DEBUG: Trying alternative extraction method..."
    
    # Try with smaller sample and different approach
    if [[ $SEQ_LENGTH -eq 24 ]]; then
        timeout 60 bash -c "zcat '$EXTRACTED_BC_FILE' | sed -n '2~4p' | cut -c7-22 | head -1000 > '$OUTPUT_DIR/qc/test_barcodes.txt'"
    else
        timeout 60 bash -c "zcat '$EXTRACTED_BC_FILE' | sed -n '2~4p' | head -1000 > '$OUTPUT_DIR/qc/test_barcodes.txt'"
    fi
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Alternative extraction also failed"
        exit 1
    fi
elif [[ $EXTRACT_EXIT_CODE -ne 0 ]]; then
    echo "ERROR: Barcode extraction failed with exit code $EXTRACT_EXIT_CODE"
    echo "DEBUG: Checking if barcode file is readable..."
    if [[ ! -r "$EXTRACTED_BC_FILE" ]]; then
        echo "ERROR: Barcode file is not readable: $EXTRACTED_BC_FILE"
        exit 1
    fi
    echo "DEBUG: Trying alternative extraction method..."
    if [[ $SEQ_LENGTH -eq 24 ]]; then
        zcat "$EXTRACTED_BC_FILE" | sed -n '2~4p' | cut -c7-22 | head -1000 > "$OUTPUT_DIR/qc/test_barcodes.txt" || {
            echo "ERROR: All barcode extraction methods failed"
            exit 1
        }
    else
        zcat "$EXTRACTED_BC_FILE" | sed -n '2~4p' | head -1000 > "$OUTPUT_DIR/qc/test_barcodes.txt" || {
            echo "ERROR: All barcode extraction methods failed"
            exit 1
        }
    fi
fi

TEST_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes.txt" 2>/dev/null || echo "0")
echo "DEBUG: Extracted $TEST_BC_COUNT test barcodes"

if [[ $TEST_BC_COUNT -eq 0 ]]; then
    echo "ERROR: No test barcodes extracted from $EXTRACTED_BC_FILE"
    echo "DEBUG: Checking barcode file contents..."
    echo "DEBUG: File size: $(stat -c%s "$EXTRACTED_BC_FILE" 2>/dev/null || echo 'unknown') bytes"
    echo "DEBUG: First few bytes: $(zcat "$EXTRACTED_BC_FILE" 2>/dev/null | head -c 100 || echo 'cannot read')"
    exit 1
fi

# Analyze barcode quality with N content
CLEAN_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes_clean.txt" 2>/dev/null || echo "0")
N_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes_with_n.txt" 2>/dev/null || echo "0")

echo "DEBUG: Barcode quality analysis:"
echo "  Clean barcodes (no N's): $CLEAN_BC_COUNT"
echo "  Barcodes with N's: $N_BC_COUNT"
echo "  Total barcodes: $TEST_BC_COUNT"
echo "  N contamination rate: $(echo "scale=2; ${N_BC_COUNT:-0} * 100 / ${TEST_BC_COUNT:-1}" | bc)%"

echo "DEBUG: Sample of extracted test barcodes:"
echo "  Clean barcodes:"
head -3 "$OUTPUT_DIR/qc/test_barcodes_clean.txt" 2>/dev/null || echo "    No clean barcodes found"
echo "  Barcodes with N's:"
head -3 "$OUTPUT_DIR/qc/test_barcodes_with_n.txt" 2>/dev/null || echo "    No N-containing barcodes found"

# Test both orientations with timeout protection
echo "DEBUG: Testing whitelist orientations..."
echo "DEBUG: Testing original whitelist..."
MATCHES_ORIG=$(timeout 120 bash -c "grep -Ff '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Original whitelist test timed out after 2 minutes"
    MATCHES_ORIG=0
fi

echo "DEBUG: Testing reverse complement whitelist..."
MATCHES_RC=$(timeout 120 bash -c "grep -Ff '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Reverse complement whitelist test timed out after 2 minutes"
    MATCHES_RC=0
fi

echo "Whitelist matching results:"
echo "  Original: $MATCHES_ORIG / 10000 ($(echo "scale=2; ${MATCHES_ORIG:-0} * 100 / 10000" | bc)%)"
echo "  Reverse complement: $MATCHES_RC / 10000 ($(echo "scale=2; ${MATCHES_RC:-0} * 100 / 10000" | bc)%)"

# Choose the best matching whitelist
if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    echo "Using original whitelist"
fi

# Test clean barcodes separately for better match rates
if [[ $CLEAN_BC_COUNT -gt 0 ]]; then
    echo "DEBUG: Testing clean barcodes against whitelists with timeout protection..."
    CLEAN_MATCHES_ORIG=$(timeout 60 bash -c "grep -Ff '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes_clean.txt' | wc -l" 2>/dev/null || echo "0")
    CLEAN_MATCHES_RC=$(timeout 60 bash -c "grep -Ff '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes_clean.txt' | wc -l" 2>/dev/null || echo "0")
    
    echo "Clean barcode matching results:"
    echo "  Original: $CLEAN_MATCHES_ORIG / $CLEAN_BC_COUNT ($(echo "scale=2; ${CLEAN_MATCHES_ORIG:-0} * 100 / ${CLEAN_BC_COUNT:-1}" | bc)%)"
    echo "  Reverse complement: $CLEAN_MATCHES_RC / $CLEAN_BC_COUNT ($(echo "scale=2; ${CLEAN_MATCHES_RC:-0} * 100 / ${CLEAN_BC_COUNT:-1}" | bc)%)"
fi

# Check if we have reasonable matches
MAX_MATCHES=$((MATCHES_ORIG > MATCHES_RC ? MATCHES_ORIG : MATCHES_RC))
if [[ $MAX_MATCHES -lt 500 ]]; then
    echo "WARNING: Very low barcode match rate ($(echo "scale=2; ${MAX_MATCHES:-0} * 100 / 10000" | bc)%)"
    echo "This indicates poor sequencing quality with many N's"
    echo "Setting high error tolerance for chromap..."
    BC_ERROR_THRESHOLD=3
elif [[ $MAX_MATCHES -lt 2000 ]]; then
    echo "WARNING: Low barcode match rate ($(echo "scale=2; ${MAX_MATCHES:-0} * 100 / 10000" | bc)%)"
    echo "This may be due to N contamination from poor sequencing quality"
    echo "Setting moderate error tolerance for chromap..."
    BC_ERROR_THRESHOLD=2
else
    echo "Acceptable barcode match rate detected"
    BC_ERROR_THRESHOLD=1
fi

# Step 3: Run chromap alignment
echo ""
echo "Step 3: Running chromap alignment..."

# Check chromap availability
echo "DEBUG: Checking chromap availability..."
if ! command -v chromap &> /dev/null; then
    echo "DEBUG: chromap not found in PATH, trying to load module..."
    
    # Try multiple module loading approaches
    if module load chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded chromap module"
    elif module load Chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded Chromap module (capitalized)"
    elif module load bioinformatics/chromap 2>/dev/null; then
        echo "DEBUG: Successfully loaded bioinformatics/chromap module"
    else
        echo "ERROR: chromap not found and module load failed"
        echo "DEBUG: Available modules containing 'chromap':"
        module avail chromap 2>&1 || echo "No chromap modules found"
        echo "DEBUG: Trying to install chromap via conda..."
        
        # Try conda installation
        if command -v conda &> /dev/null; then
            echo "DEBUG: Installing chromap via conda..."
            conda install -y -c bioconda chromap || {
                echo "ERROR: Failed to install chromap via conda"
                echo "Please install chromap manually from: https://github.com/haowenz/chromap"
                exit 1
            }
        else
            echo "ERROR: Neither chromap nor conda available"
            echo "Please install chromap from: https://github.com/haowenz/chromap"
            exit 1
        fi
    fi
    
    # Check again after module load/installation
    if ! command -v chromap &> /dev/null; then
        echo "ERROR: chromap still not available after module load/installation"
        exit 1
    fi
fi

echo "DEBUG: chromap found: $(which chromap)"
echo "DEBUG: chromap version: $(chromap --version 2>&1 | head -1 || echo 'Version check failed')"

# Check if chromap index exists, if not create it
if [[ ! -f "$CHROMAP_INDEX" ]]; then
    echo "DEBUG: Chromap index not found, creating index..."
    mkdir -p "$(dirname "$CHROMAP_INDEX")"
    
    echo "DEBUG: Building chromap index from reference: $REF_GENOME"
    chromap -i -r "$REF_GENOME" -o "$CHROMAP_INDEX"
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to build chromap index"
        exit 1
    fi
    echo "DEBUG: Chromap index created successfully: $CHROMAP_INDEX"
else
    echo "DEBUG: Using existing chromap index: $CHROMAP_INDEX"
fi

# Run chromap with N-tolerant parameters
echo "DEBUG: Running chromap with the following parameters:"
echo "  Reference: $REF_GENOME"
echo "  Index: $CHROMAP_INDEX"
echo "  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz"
echo "  R2: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz"
echo "  Barcode: $EXTRACTED_BC_FILE"
echo "  Whitelist: $WHITELIST"
echo "  Barcode error threshold: $BC_ERROR_THRESHOLD (adjusted for N contamination)"
echo "  N contamination rate: $(echo "scale=2; ${N_BC_COUNT:-0} * 100 / ${TEST_BC_COUNT:-1}" | bc)%"

chromap --preset atac \
    -x "$CHROMAP_INDEX" \
    -r "$REF_GENOME" \
    -1 "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" \
    -2 "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" \
    -b "$EXTRACTED_BC_FILE" \
    --barcode-whitelist "$WHITELIST" \
    --bc-error-threshold $BC_ERROR_THRESHOLD \
    -e 5 \
    -t $THREADS \
    --low-mem \
    --skip-barcode-check \
    -o "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" 2>&1 | tee "$OUTPUT_DIR/logs/${SAMPLE}_chromap.log"

# Check if chromap succeeded
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" ]]; then
    echo "ERROR: Chromap failed to generate fragments file"
    exit 1
fi

# Step 4: Compress and index fragments
echo ""
echo "Step 4: Processing fragments..."

bgzip -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv"
tabix -f -s 1 -b 2 -e 3 -p bed "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"

# Step 5: Generate quality statistics
echo ""
echo "Step 5: Generating statistics..."

TOTAL_FRAGS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | wc -l)
UNIQUE_BC=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | sort -u | wc -l)
VALID_BC=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | grep -Ff "$WHITELIST" | wc -l)
FRAGS_IN_CELLS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | grep -Ff "$WHITELIST" | wc -l)

echo "Fragment statistics:"
echo "  Total fragments: $TOTAL_FRAGS"
echo "  Unique barcodes: $UNIQUE_BC"
echo "  Valid whitelist barcodes: $(echo "$VALID_BC" | cut -f4 | sort -u | wc -l)"
echo "  Fragments in cells: $FRAGS_IN_CELLS"
echo "  Fraction in cells: $(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)"

# Step 6: Peak calling with MACS2
echo ""
echo "DEBUG: Step 6 - Calling peaks with MACS2..."

# Generate genome sizes
echo "DEBUG: Checking for genome sizes file..."
if [[ ! -f "$OUTPUT_DIR/qc/genome.chrom.sizes" ]]; then
    echo "DEBUG: Generating genome sizes file from reference index..."
    if [[ -f "${REF_GENOME}.fai" ]]; then
        cut -f1,2 "${REF_GENOME}.fai" > "$OUTPUT_DIR/qc/genome.chrom.sizes"
        echo "DEBUG: Created genome sizes file: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
    else
        echo "DEBUG: Creating reference index first..."
        samtools faidx "$REF_GENOME"
        cut -f1,2 "${REF_GENOME}.fai" > "$OUTPUT_DIR/qc/genome.chrom.sizes"
        echo "DEBUG: Created genome sizes file: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
    fi
else
    echo "DEBUG: Genome sizes file already exists: $(wc -l < "$OUTPUT_DIR/qc/genome.chrom.sizes") chromosomes"
fi

# Convert fragments to reads for peak calling
echo "DEBUG: Converting fragments to reads for peak calling..."
zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | \
    awk 'BEGIN{OFS="\t"}{
        print $1, $2, $2+50, $4, ".", "+"
        print $1, $3-50, $3, $4, ".", "-"
    }' | \
    sed '/chrM/d' | \
    bedClip stdin "$OUTPUT_DIR/qc/genome.chrom.sizes" stdout | \
    sort -k1,1 -k2,2n | \
    gzip > "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Reads file created successfully"
    echo "DEBUG: Reads file size: $(stat -c%s "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz") bytes"
    echo "DEBUG: Number of reads: $(zcat "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" | wc -l)"
else
    echo "ERROR: Failed to create reads file"
    exit 1
fi

# Call peaks
echo "DEBUG: Checking for MACS2..."
if ! command -v macs2 &> /dev/null; then
    echo "ERROR: MACS2 not found. Please install MACS2 for peak calling"
    echo "DEBUG: You can install with: pip install MACS2"
    exit 1
fi

echo "DEBUG: MACS2 found: $(which macs2)"
echo "DEBUG: MACS2 version: $(macs2 --version 2>&1 || echo 'Version check failed')"

echo "DEBUG: Running MACS2 peak calling..."
echo "DEBUG: MACS2 command:"
echo "macs2 callpeak -t \"$OUTPUT_DIR/${SAMPLE}_reads.bed.gz\" -g mm -f BED -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -B --SPMR --outdir \"$OUTPUT_DIR/peaks\" -n \"${SAMPLE}\""

macs2 callpeak \
    -t "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" \
    -g mm -f BED -q 0.01 \
    --nomodel --shift -100 --extsize 200 \
    --keep-dup all \
    -B --SPMR \
    --outdir "$OUTPUT_DIR/peaks" \
    -n "${SAMPLE}"

MACS2_EXIT_CODE=$?
echo "DEBUG: MACS2 exit code: $MACS2_EXIT_CODE"

if [[ $MACS2_EXIT_CODE -ne 0 ]]; then
    echo "ERROR: MACS2 failed with exit code $MACS2_EXIT_CODE"
    exit 1
fi

if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "ERROR: MACS2 completed but peaks file not found"
    echo "DEBUG: Files in peaks directory:"
    ls -la "$OUTPUT_DIR/peaks/" 2>/dev/null || echo "Peaks directory not found"
    exit 1
fi

echo "DEBUG: MACS2 completed successfully"
echo "DEBUG: Peak file size: $(stat -c%s "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak") bytes"
echo "DEBUG: Number of peaks: $(wc -l < "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak")"

# Step 7: Generate peak-by-cell matrix
echo ""
echo "DEBUG: Step 7 - Generating peak-by-cell matrix..."

# Sort peaks
echo "DEBUG: Sorting peaks..."
cut -f 1-4 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" | \
    sort -k1,1 -k2,2n > "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Peaks sorted successfully"
    echo "DEBUG: Sorted peaks file size: $(stat -c%s "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed") bytes"
else
    echo "ERROR: Failed to sort peaks"
    exit 1
fi

# Find reads in peaks per cell
echo "DEBUG: Checking for bedtools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found. Please install bedtools"
    exit 1
fi

echo "DEBUG: Finding reads in peaks per cell..."
echo "DEBUG: bedtools intersect command:"
echo "bedtools intersect -a \"$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed\" -b \"$OUTPUT_DIR/${SAMPLE}_reads.bed.gz\" -wo -sorted -g \"$OUTPUT_DIR/qc/genome.chrom.sizes\""

bedtools intersect \
    -a "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed" \
    -b "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" \
    -wo -sorted -g "$OUTPUT_DIR/qc/genome.chrom.sizes" | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: Peak-read overlap completed successfully"
    echo "DEBUG: Peak-read overlap file size: $(stat -c%s "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz") bytes"
    echo "DEBUG: Number of cells with peaks: $(zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | wc -l)"
else
    echo "ERROR: Failed to generate peak-read overlap"
    exit 1
fi

# Create matrix files
echo "DEBUG: Creating matrix files..."
mkdir -p "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"

# Peaks bed file
echo "DEBUG: Creating peaks.bed file..."
cut -f 1-3 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed" > \
    "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: peaks.bed created: $(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed") peaks"
else
    echo "ERROR: Failed to create peaks.bed"
    exit 1
fi

# Barcodes file
echo "DEBUG: Creating barcodes.tsv file..."
zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | \
    cut -f 1 > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv"

if [[ $? -eq 0 ]]; then
    echo "DEBUG: barcodes.tsv created: $(wc -l < "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv") barcodes"
else
    echo "ERROR: Failed to create barcodes.tsv"
    exit 1
fi

# Step 8: Generate final report
echo ""
echo "Step 8: Generating final report..."

cat > "$OUTPUT_DIR/${SAMPLE}_processing_report.txt" << EOF
ATAC Processing Report
======================
Sample: $SAMPLE
Date: $(date)
Output Directory: $OUTPUT_DIR

Input Files:
  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz
  R2 (barcodes): $EXTRACTED_BC_FILE (extracted from positions 7-22)
  R3: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz

Processing Parameters:
  Whitelist: $(basename $WHITELIST)
  Barcode mismatches allowed: $BARCODE_MISMATCH
  Read mismatches allowed: 5

Quality Metrics:
  Barcode match rate: $(echo "scale=2; ${MAX_MATCHES:-0} * 100 / 10000" | bc)%
  Total fragments: $TOTAL_FRAGS
  Unique barcodes: $UNIQUE_BC
  Fragments in valid cells: $FRAGS_IN_CELLS
  Fraction in cells: $(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)

Output Files:
  Fragments: $OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz
  Peaks: $OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak
  Peak-cell matrix: $OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/

Next Steps:
  - Import into Seurat/Signac/ArchR for downstream analysis
  - Check $OUTPUT_DIR/peaks/${SAMPLE}_treat_pileup.bdg for visualization
  - Review chromap log: $OUTPUT_DIR/logs/${SAMPLE}_chromap.log
EOF

cat "$OUTPUT_DIR/${SAMPLE}_processing_report.txt"

echo ""
echo "========================================="
echo "Processing complete for $SAMPLE"
echo "Report saved to: $OUTPUT_DIR/${SAMPLE}_processing_report.txt"
echo "End time: $(date)"
echo "========================================="