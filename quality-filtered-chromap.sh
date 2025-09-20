#!/bin/bash
#SBATCH --job-name=atac_quality_filter
#SBATCH --output=logs/atac_quality_%A.out
#SBATCH --error=logs/atac_quality_%A.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

# Quality-filtered ATAC processing pipeline
# Handles barcodes with N bases and poor quality sequences

set -euo pipefail

# Configuration
SAMPLE="R26-Nestin-Mut-adult"  # Modify as needed
DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/quality_filtered_output"
REF_GENOME="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/mm10/mm10.fa"
CHROMAP_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/mm10/chromap_index/genome.index"
THREADS=16

mkdir -p "$OUTPUT_DIR/logs" "$OUTPUT_DIR/qc"

echo "========================================="
echo "Processing sample: $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# Step 1: Assess data quality
echo "Step 1: Assessing data quality..."

# Check R2 quality (barcodes)
echo "Checking barcode quality in R2..."
zcat "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" | \
    awk 'NR%4==2 {print substr($0,1,16)}' | \
    head -100000 > "$OUTPUT_DIR/qc/barcodes_100k.txt"

# Count N bases
TOTAL_BC=$(wc -l < "$OUTPUT_DIR/qc/barcodes_100k.txt")
N_BC=$(grep -c 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" || true)
CLEAN_BC=$(grep -vc 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" || true)
N_PERCENT=$(echo "scale=2; $N_BC * 100 / $TOTAL_BC" | bc)

echo "Barcode quality summary (first 100k reads):"
echo "  Total barcodes: $TOTAL_BC"
echo "  Barcodes with N's: $N_BC ($N_PERCENT%)"
echo "  Clean barcodes: $CLEAN_BC"

# Get whitelist files
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/barcode_whitelists"
mkdir -p "$WHITELIST_DIR"

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1.txt" ]]; then
    echo "Downloading barcode whitelist..."
    wget -q -O "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz" \
        https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
    gunzip "$WHITELIST_DIR/atac_737K-arc-v1.txt.gz"
fi

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

# Check orientation using only clean barcodes
grep -v 'N' "$OUTPUT_DIR/qc/barcodes_100k.txt" | head -10000 > "$OUTPUT_DIR/qc/clean_barcodes_10k.txt"
MATCHES_ORIG=$(grep -Ff "$WHITELIST_DIR/atac_737K-arc-v1.txt" "$OUTPUT_DIR/qc/clean_barcodes_10k.txt" | wc -l)
MATCHES_RC=$(grep -Ff "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" "$OUTPUT_DIR/qc/clean_barcodes_10k.txt" | wc -l)

echo "Barcode whitelist matching (clean barcodes only):"
echo "  Original whitelist: $MATCHES_ORIG / 10000"
echo "  Reverse complement: $MATCHES_RC / 10000"

if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    echo "Using original whitelist"
fi

# Step 2: Quality filtering with fastp
echo ""
echo "Step 2: Quality filtering with fastp..."

# Check if fastp is available
if ! command -v fastp &> /dev/null; then
    echo "fastp not found. Installing via conda or loading module..."
    # Try to load module first
    module load fastp 2>/dev/null || true
    
    # If still not available, install with conda
    if ! command -v fastp &> /dev/null; then
        echo "Installing fastp..."
        conda install -y -c bioconda fastp || {
            echo "ERROR: Cannot install fastp. Please install manually."
            exit 1
        }
    fi
fi

# Filter all reads with fastp
# This will remove reads with too many N's and low quality
echo "Filtering R1 (genomic read 1)..."
fastp -i "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" \
      -o "$OUTPUT_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
      -h "$OUTPUT_DIR/qc/${SAMPLE}_R1_fastp.html" \
      -j "$OUTPUT_DIR/qc/${SAMPLE}_R1_fastp.json" \
      --thread $THREADS \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --n_base_limit 5 \
      --length_required 36

echo "Filtering R2 (barcode)..."
# For barcodes, be less strict since we need the first 16bp
fastp -i "$DATA_DIR/${SAMPLE}_R2_001.fastq.gz" \
      -o "$OUTPUT_DIR/${SAMPLE}_R2_filtered.fastq.gz" \
      -h "$OUTPUT_DIR/qc/${SAMPLE}_R2_fastp.html" \
      -j "$OUTPUT_DIR/qc/${SAMPLE}_R2_fastp.json" \
      --thread $THREADS \
      --qualified_quality_phred 15 \
      --unqualified_percent_limit 50 \
      --n_base_limit 2 \
      --length_required 24 \
      --max_len1 24

echo "Filtering R3 (genomic read 2)..."
fastp -i "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" \
      -o "$OUTPUT_DIR/${SAMPLE}_R3_filtered.fastq.gz" \
      -h "$OUTPUT_DIR/qc/${SAMPLE}_R3_fastp.html" \
      -j "$OUTPUT_DIR/qc/${SAMPLE}_R3_fastp.json" \
      --thread $THREADS \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --n_base_limit 5 \
      --length_required 36

# Step 3: Extract 16bp barcodes
echo ""
echo "Step 3: Extracting 16bp barcodes..."

zcat "$OUTPUT_DIR/${SAMPLE}_R2_filtered.fastq.gz" | \
    awk 'NR%4==1 {print $0} 
         NR%4==2 {bc=substr($0,1,16); if(bc !~ /N/) print bc; else print "NNNNNNNNNNNNNNNN"} 
         NR%4==3 {print "+"} 
         NR%4==0 {print substr($0,1,16)}' | \
    gzip > "$OUTPUT_DIR/${SAMPLE}_R2_16bp.fastq.gz"

# Step 4: Create barcode error correction map
echo ""
echo "Step 4: Creating barcode error correction map..."

# Extract all barcodes that have 0 or 1 N
zcat "$OUTPUT_DIR/${SAMPLE}_R2_16bp.fastq.gz" | \
    awk 'NR%4==2' | \
    awk '{n=gsub(/N/,"N"); if(n<=1) print}' | \
    sort | uniq -c | sort -rn | \
    awk '{print $2}' > "$OUTPUT_DIR/qc/observed_barcodes.txt"

# For barcodes with 1 N, try to correct to nearest whitelist barcode
echo "Attempting barcode error correction..."
python3 - << 'EOF' "$WHITELIST" "$OUTPUT_DIR/qc/observed_barcodes.txt" "$OUTPUT_DIR/barcode_correction.txt"
import sys
from collections import defaultdict

whitelist_file = sys.argv[1]
observed_file = sys.argv[2]
output_file = sys.argv[3]

# Load whitelist
whitelist = set()
with open(whitelist_file) as f:
    for line in f:
        whitelist.add(line.strip())

# Process observed barcodes
corrections = {}
with open(observed_file) as f:
    for line in f:
        bc = line.strip()
        if 'N' not in bc:
            if bc in whitelist:
                corrections[bc] = bc  # Valid barcode
        elif bc.count('N') == 1:
            # Try all 4 possible substitutions
            n_pos = bc.index('N')
            for base in 'ACGT':
                candidate = bc[:n_pos] + base + bc[n_pos+1:]
                if candidate in whitelist:
                    corrections[bc] = candidate
                    break

# Write corrections
with open(output_file, 'w') as f:
    for orig, corrected in corrections.items():
        f.write(f"{orig}\t{corrected}\n")

print(f"Created correction map for {len(corrections)} barcodes")
EOF

# Step 5: Apply barcode corrections and filter
echo ""
echo "Step 5: Applying barcode corrections..."

# This step would require custom code to replace barcodes in the FASTQ
# For now, we'll proceed with chromap which can handle some mismatches

# Step 6: Run chromap with filtered data
echo ""
echo "Step 6: Running chromap alignment..."

# Check chromap availability
if ! command -v chromap &> /dev/null; then
    echo "ERROR: chromap not found. Please install chromap first."
    echo "You can download from: https://github.com/haowenz/chromap"
    exit 1
fi

chromap --preset atac \
    -x "$CHROMAP_INDEX" \
    -r "$REF_GENOME" \
    -1 "$OUTPUT_DIR/${SAMPLE}_R1_filtered.fastq.gz" \
    -2 "$OUTPUT_DIR/${SAMPLE}_R3_filtered.fastq.gz" \
    -b "$OUTPUT_DIR/${SAMPLE}_R2_16bp.fastq.gz" \
    --barcode-whitelist "$WHITELIST" \
    --barcode-mismatch 1 \
    --read-mismatch 4 \
    --threads $THREADS \
    -o "$OUTPUT_DIR/${SAMPLE}_fragments.tsv" \
    --low-mem

# Step 7: Post-processing
echo ""
echo "Step 7: Post-processing..."

# Compress and index
bgzip -f "$OUTPUT_DIR/${SAMPLE}_fragments.tsv"
tabix -f -s 1 -b 2 -e 3 -p bed "$OUTPUT_DIR/${SAMPLE}_fragments.tsv.gz"

# Count statistics
TOTAL_FRAGS=$(zcat "$OUTPUT_DIR/${SAMPLE}_fragments.tsv.gz" | wc -l)
UNIQUE_BC=$(zcat "$OUTPUT_DIR/${SAMPLE}_fragments.tsv.gz" | cut -f4 | sort -u | wc -l)
VALID_BC=$(zcat "$OUTPUT_DIR/${SAMPLE}_fragments.tsv.gz" | cut -f4 | grep -Ff "$WHITELIST" | sort -u | wc -l)

echo ""
echo "========================================="
echo "Final Statistics:"
echo "  Total fragments: $TOTAL_FRAGS"
echo "  Unique barcodes: $UNIQUE_BC"
echo "  Valid whitelist barcodes: $VALID_BC"
echo "  Valid barcode rate: $(echo "scale=2; $VALID_BC * 100 / $UNIQUE_BC" | bc)%"
echo "========================================="

# Step 8: Generate quality report
echo ""
echo "Step 8: Generating quality report..."

cat > "$OUTPUT_DIR/quality_report_${SAMPLE}.txt" << EOL
ATAC Processing Quality Report
Sample: $SAMPLE
Date: $(date)

Input Data Quality:
  Barcodes with N bases: $N_PERCENT%
  Clean barcode rate: $(echo "scale=2; $CLEAN_BC * 100 / $TOTAL_BC" | bc)%

Whitelist Matching:
  Whitelist used: $(basename "$WHITELIST")
  Match rate (clean barcodes): $(echo "scale=2; $MATCHES_ORIG * 100 / 10000" | bc)%

Output Statistics:
  Total fragments: $TOTAL_FRAGS
  Unique barcodes: $UNIQUE_BC
  Valid barcodes: $VALID_BC
  Valid rate: $(echo "scale=2; $VALID_BC * 100 / $UNIQUE_BC" | bc)%

Quality Filtering:
  See HTML reports in qc/ directory for detailed statistics

Recommendations:
$(if [[ $(echo "$N_PERCENT > 20" | bc) -eq 1 ]]; then
    echo "  - High N base rate indicates sequencing quality issues"
    echo "  - Consider re-sequencing or more aggressive filtering"
fi)
$(if [[ $VALID_BC -lt 1000 ]]; then
    echo "  - Low number of valid barcodes detected"
    echo "  - Check if correct chemistry/whitelist is being used"
fi)
EOL

echo ""
echo "Processing complete! Check $OUTPUT_DIR/quality_report_${SAMPLE}.txt for summary"
echo "End time: $(date)"