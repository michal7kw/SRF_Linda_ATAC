#!/bin/bash
#SBATCH --job-name=atac_processing_fixed
#SBATCH --output=logs/atac_fixed_%a.out
#SBATCH --error=logs/atac_fixed_%a.err
#SBATCH --array=0
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

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
echo "Processing ATAC sample: $SAMPLE (FIXED VERSION)"
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

# Step 1: Extract barcodes from RIGHTMOST 16bp of 24bp R2 (positions 9-24)
echo "Step 1: Extracting barcodes from rightmost 16bp of R2..."
EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

echo "DEBUG: R2 structure according to facility:"
echo "  [8bp SPACER][16bp 10x barcode] = 24bp total"
echo "  Extracting rightmost 16bp (positions 9-24)"

# Ensure any previous partial file is removed to prevent errors
# rm -f "$EXTRACTED_BC_FILE"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "DEBUG: Extracting rightmost 16bp from R2..."
    mkdir -p "$PROCESSED_DATA_DIR"
    
    # Extract rightmost 16bp from 24bp R2 read
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,9,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,9,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
        
    echo "DEBUG: Barcode extraction completed"
else
    echo "DEBUG: Using existing extracted barcodes: $EXTRACTED_BC_FILE"
fi

echo "Using extracted barcodes from: $EXTRACTED_BC_FILE"

# Step 2: Test barcode extraction quality
echo "Step 2: Testing barcode extraction quality..."
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
mkdir -p "$WHITELIST_DIR"

# Download whitelists if needed
if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1.txt" ]]; then
    echo "DEBUG: Downloading whitelist..."
    wget -q -O "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to download whitelist"
        exit 1
    fi
    
    gunzip "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
    echo "DEBUG: Whitelist downloaded and extracted"
fi

echo "DEBUG: Checking whitelist files..."
echo "  Original whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1.txt") barcodes"

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    echo "DEBUG: Creating reverse complement whitelist..."
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

echo "  RC whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt") barcodes"

# Test extracted barcodes
echo "DEBUG: Testing extracted barcodes against whitelists..."
echo "DEBUG: Extracting test barcodes..."

# Extract test barcodes (should now be exactly 16bp)
timeout 300 bash -c "zcat '$EXTRACTED_BC_FILE' | awk 'NR%4==2' | head -10000 > '$OUTPUT_DIR/qc/test_barcodes.txt'"
EXTRACT_EXIT_CODE=$?

if [[ $EXTRACT_EXIT_CODE -eq 124 ]]; then
    echo "ERROR: Barcode extraction timed out"
    exit 1
elif [[ $EXTRACT_EXIT_CODE -ne 0 ]]; then
    echo "ERROR: Barcode extraction failed"
    exit 1
fi

TEST_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes.txt")
echo "DEBUG: Extracted $TEST_BC_COUNT test barcodes"

if [[ $TEST_BC_COUNT -eq 0 ]]; then
    echo "ERROR: No test barcodes extracted"
    exit 1
fi

# Check barcode length
SAMPLE_SEQ=$(head -1 "$OUTPUT_DIR/qc/test_barcodes.txt")
SEQ_LENGTH=${#SAMPLE_SEQ}
echo "DEBUG: Sample barcode: $SAMPLE_SEQ"
echo "DEBUG: Barcode length: $SEQ_LENGTH"

if [[ $SEQ_LENGTH -ne 16 ]]; then
    echo "ERROR: Expected 16bp barcodes, got ${SEQ_LENGTH}bp"
    echo "DEBUG: This suggests the R2 structure may not be as expected"
    exit 1
fi

# Analyze N content
N_COUNT=$(echo "$SAMPLE_SEQ" | grep -o 'N' | wc -l)
echo "DEBUG: Number of N's in sample barcode: $N_COUNT"

# Test both whitelist orientations with timeout and optimization
echo "DEBUG: Testing whitelist orientations..."
echo "DEBUG: Using optimized matching approach to avoid hangs..."

# Use a smaller test set for initial orientation testing
head -1000 "$OUTPUT_DIR/qc/test_barcodes.txt" > "$OUTPUT_DIR/qc/test_barcodes_small.txt"

echo "DEBUG: Testing original whitelist..."
MATCHES_ORIG=$(timeout 60 bash -c "grep -Fxf '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Original whitelist test timed out, using alternative method"
    MATCHES_ORIG=$(timeout 30 bash -c "awk 'NR==FNR{a[\$0]=1;next} \$0 in a{c++} END{print c+0}' '$WHITELIST_DIR/atac_737K-arc-v1.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt'" 2>/dev/null || echo "0")
fi

echo "DEBUG: Testing reverse complement whitelist..."
MATCHES_RC=$(timeout 60 bash -c "grep -Fxf '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt' | wc -l" 2>/dev/null || echo "0")
if [[ $? -eq 124 ]]; then
    echo "WARNING: Reverse complement whitelist test timed out, using alternative method"
    MATCHES_RC=$(timeout 30 bash -c "awk 'NR==FNR{a[\$0]=1;next} \$0 in a{c++} END{print c+0}' '$WHITELIST_DIR/atac_737K-arc-v1_rc.txt' '$OUTPUT_DIR/qc/test_barcodes_small.txt'" 2>/dev/null || echo "0")
fi

SMALL_BC_COUNT=$(wc -l < "$OUTPUT_DIR/qc/test_barcodes_small.txt")
echo "Whitelist matching results (tested on $SMALL_BC_COUNT barcodes):"
echo "  Original: $MATCHES_ORIG / $SMALL_BC_COUNT ($(echo "scale=2; ${MATCHES_ORIG:-0} * 100 / $SMALL_BC_COUNT" | bc)%)"
echo "  Reverse complement: $MATCHES_RC / $SMALL_BC_COUNT ($(echo "scale=2; ${MATCHES_RC:-0} * 100 / $SMALL_BC_COUNT" | bc)%)"

# Choose the best matching whitelist
if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    echo "Using original whitelist"
fi

# Set barcode error threshold based on match rate
MAX_MATCHES=$((MATCHES_ORIG > MATCHES_RC ? MATCHES_ORIG : MATCHES_RC))
MATCH_RATE=$(echo "scale=2; ${MAX_MATCHES:-0} * 100 / $SMALL_BC_COUNT" | bc)

echo "DEBUG: Best match rate: ${MATCH_RATE}% (from $SMALL_BC_COUNT test barcodes)"

if [[ $(echo "$MATCH_RATE < 5.0" | bc) -eq 1 ]]; then
    echo "ERROR: Very low barcode match rate (${MATCH_RATE}%)"
    echo "This suggests either:"
    echo "  1. Wrong extraction positions (check R2 structure)"
    echo "  2. Wrong whitelist file"
    echo "  3. Severely degraded sequencing quality"
    exit 1
elif [[ $(echo "$MATCH_RATE < 20.0" | bc) -eq 1 ]]; then
    echo "WARNING: Low barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=2
elif [[ $(echo "$MATCH_RATE < 50.0" | bc) -eq 1 ]]; then
    echo "INFO: Moderate barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=1
else
    echo "INFO: Good barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=1
fi

# Step 2.5: Validate read counts match between files
echo ""
echo "Step 2.5: Validating read counts across files..."

echo "DEBUG: Counting reads in each file..."
R1_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" | wc -l | awk '{print $1/4}')
R3_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" | wc -l | awk '{print $1/4}') 
BC_COUNT=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')

echo "Read counts:"
echo "  R1 (genomic, forward): $R1_COUNT"
echo "  R3 (genomic, reverse): $R3_COUNT" 
echo "  R2 (extracted barcodes): $BC_COUNT"

# Check if all counts match
if [[ "$R1_COUNT" -ne "$R3_COUNT" ]]; then
    echo "ERROR: R1 and R3 read counts don't match ($R1_COUNT vs $R3_COUNT)"
    exit 1
fi

if [[ "$R1_COUNT" -ne "$BC_COUNT" ]]; then
    echo "ERROR: R1 and barcode counts don't match ($R1_COUNT vs $BC_COUNT)"
    echo "This suggests the R2 barcode file was truncated or has different structure"
    echo "Attempting to fix by re-extracting barcodes..."
    
    # Re-extract barcodes ensuring same read count as R1
    echo "DEBUG: Re-extracting barcodes with read count validation..."
    ORIG_R2_COUNT=$(zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | wc -l | awk '{print $1/4}')
    echo "DEBUG: Original R2 read count: $ORIG_R2_COUNT"
    
    if [[ "$ORIG_R2_COUNT" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Original R2 file has different read count than R1 ($ORIG_R2_COUNT vs $R1_COUNT)"
        echo "This indicates corrupted or incomplete FASTQ files"
        exit 1
    fi
    
    # Re-extract with validation
    rm -f "$EXTRACTED_BC_FILE"
    zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
        awk '{
            if(NR%4==1) print $0;
            else if(NR%4==2) print substr($0,9,16);
            else if(NR%4==3) print $0;
            else if(NR%4==0) print substr($0,9,16);
        }' | gzip > "$EXTRACTED_BC_FILE"
    
    # Recount after extraction
    BC_COUNT_NEW=$(zcat "$EXTRACTED_BC_FILE" | wc -l | awk '{print $1/4}')
    echo "DEBUG: New barcode count after re-extraction: $BC_COUNT_NEW"
    
    if [[ "$BC_COUNT_NEW" -ne "$R1_COUNT" ]]; then
        echo "ERROR: Re-extraction failed, counts still don't match ($BC_COUNT_NEW vs $R1_COUNT)"
        exit 1
    fi
    
    echo "DEBUG: Re-extraction successful, all counts now match"
    BC_COUNT=$BC_COUNT_NEW
fi

echo "✓ All file read counts match: $R1_COUNT reads per file"

# Step 3: Run chromap alignment
echo ""
echo "Step 3: Running chromap alignment..."
echo "DEBUG: Proceeding with barcode error threshold: $BC_ERROR_THRESHOLD"

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

# Run chromap with corrected barcode file
echo "DEBUG: Running chromap with the following parameters:"
echo "  Reference: $REF_GENOME"
echo "  Index: $CHROMAP_INDEX"
echo "  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz"
echo "  R2: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz"
echo "  Barcode: $EXTRACTED_BC_FILE (rightmost 16bp from R2)"
echo "  Whitelist: $WHITELIST"
echo "  Barcode error threshold: $BC_ERROR_THRESHOLD"
echo "  Expected match rate: ${MATCH_RATE}%"

if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" ]]; then
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
        -o "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" 2>&1 | tee "$OUTPUT_DIR/logs/${SAMPLE}_chromap.log"
else
    echo "DEBUG: Skipping chromap alignment, fragments file already exists."
fi

# Check if chromap succeeded
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv" && ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" ]]; then
    echo "ERROR: Chromap failed to generate fragments file"
    echo "DEBUG: Check chromap log: $OUTPUT_DIR/logs/${SAMPLE}_chromap.log"
    exit 1
fi

# Step 4: Compress and index fragments
echo ""
echo "Step 4: Processing fragments..."

if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz.tbi" ]]; then
    bgzip -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv"
    tabix -f -s 1 -b 2 -e 3 -p bed "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
else
    echo "DEBUG: Skipping fragment indexing, index file already exists."
fi

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

if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    # Check if bedClip is available
    if command -v bedClip &> /dev/null; then
        echo "DEBUG: Using bedClip for coordinate clipping"
        zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | \
            awk 'BEGIN{OFS="\t"}{
                print $1, $2, $2+50, $4, ".", "+"
                print $1, $3-50, $3, $4, ".", "-"
            }' | \
            sed '/chrM/d' | \
            bedClip stdin "$OUTPUT_DIR/qc/genome.chrom.sizes" stdout | \
            sort -k1,1 -k2,2n | \
            gzip > "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
    else
        echo "DEBUG: bedClip not found, using alternative coordinate clipping"
        # Alternative method using awk to clip coordinates
        zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | \
            awk -v sizefile="$OUTPUT_DIR/qc/genome.chrom.sizes" '
            BEGIN {
                OFS="\t"
                # Load chromosome sizes
                while((getline line < sizefile) > 0) {
                    split(line, a, "\t")
                    chrsize[a[1]] = a[2]
                }
                close(sizefile)
            }
            {
                # Skip mitochondrial reads
                if($1 == "chrM") next
                
                # Generate reads from fragment ends
                start1 = $2; end1 = $2 + 50
                start2 = $3 - 50; end2 = $3
                
                # Clip to chromosome boundaries
                if(start1 < 0) start1 = 0
                if(end1 > chrsize[$1]) end1 = chrsize[$1]
                if(start2 < 0) start2 = 0
                if(end2 > chrsize[$1]) end2 = chrsize[$1]
                
                # Only output if coordinates are valid
                if(start1 < end1) print $1, start1, end1, $4, ".", "+"
                if(start2 < end2) print $1, start2, end2, $4, ".", "-"
            }' | \
            sort -k1,1 -k2,2n | \
            gzip > "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
    fi

    if [[ $? -eq 0 ]]; then
        echo "DEBUG: Reads file created successfully"
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

    echo "DEBUG: Running MACS2 peak calling..."
    macs2 callpeak \
        -t "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" \
        -g mm -f BED -q 0.01 \
        --nomodel --shift -100 --extsize 200 \
        --keep-dup all \
        -B --SPMR \
        --outdir "$OUTPUT_DIR/peaks" \
        -n "${SAMPLE}"
else
    echo "DEBUG: Skipping MACS2 peak calling, peaks file already exists."
fi

MACS2_EXIT_CODE=$?
if [[ $MACS2_EXIT_CODE -ne 0 ]]; then
    echo "ERROR: MACS2 failed with exit code $MACS2_EXIT_CODE"
    exit 1
fi

if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "ERROR: MACS2 completed but peaks file not found"
    exit 1
fi

echo "DEBUG: Number of peaks: $(wc -l < "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak")"

# Step 7: Generate peak-by-cell matrix
echo ""
echo "DEBUG: Step 7 - Generating peak-by-cell matrix..."

# Sort peaks
if [[ ! -f "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" ]]; then
    cut -f 1-4 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" | \
        sort -k1,1 -k2,2n > "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

    # Find reads in peaks per cell
    echo "DEBUG: Checking for bedtools..."
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: bedtools not found. Please install bedtools"
        exit 1
    fi

    bedtools intersect \
        -a "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed" \
        -b "$OUTPUT_DIR/${SAMPLE}_reads.bed.gz" \
        -wo -sorted -g "$OUTPUT_DIR/qc/genome.chrom.sizes" | \
        sort -k8,8 | \
        bedtools groupby -g 8 -c 4 -o freqdesc | \
        gzip > "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz"
else
    echo "DEBUG: Skipping peak-by-cell matrix generation, file already exists."
fi

# Create matrix files
mkdir -p "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"

# Peaks bed file
cut -f 1-3 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed" > \
    "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/peaks.bed"

# Barcodes file
zcat "$OUTPUT_DIR/${SAMPLE}_peak_read_ov.tsv.gz" | \
    cut -f 1 > "$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix/barcodes.tsv"

# Step 8: Generate final report
echo ""
echo "Step 8: Generating final report..."

cat > "$OUTPUT_DIR/${SAMPLE}_processing_report_fixed.txt" << EOF
ATAC Processing Report (FIXED VERSION)
======================================
Sample: $SAMPLE
Date: $(date)
Output Directory: $OUTPUT_DIR

CORRECTION APPLIED:
  Previous extraction: positions 7-22 from R2 (incorrect)
  Fixed extraction: rightmost 16bp from R2 (positions 9-24)
  R2 structure: [8bp SPACER][16bp 10x barcode] = 24bp total

Input Files:
  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz
  R2 (barcodes): $EXTRACTED_BC_FILE (rightmost 16bp from R2)
  R3: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz

Processing Parameters:
  Whitelist: $(basename $WHITELIST)
  Barcode error threshold: $BC_ERROR_THRESHOLD
  Read mismatches allowed: 5

Quality Metrics:
  Barcode match rate: ${MATCH_RATE}%
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

cat "$OUTPUT_DIR/${SAMPLE}_processing_report_fixed.txt"

echo ""
echo "========================================="
echo "Processing complete for $SAMPLE (FIXED VERSION)"
echo "Report saved to: $OUTPUT_DIR/${SAMPLE}_processing_report_fixed.txt"
echo "End time: $(date)"
echo "========================================="