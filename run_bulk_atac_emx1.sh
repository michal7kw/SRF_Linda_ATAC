#!/bin/bash
#SBATCH --job-name=bulk_atac_emx1
#SBATCH --output=logs/bulk_atac_emx1_%A_%a.out
#SBATCH --error=logs/bulk_atac_emx1_%A_%a.err
#SBATCH --array=0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=workq

# Bulk ATAC-seq Processing Pipeline for Emx1 Sample
# This script processes ATAC-seq data as bulk, bypassing single-cell barcode issues.
# It performs trimming, alignment, filtering, and peak calling.

set -euo pipefail

# Cleanup function
cleanup() {
    EXIT_CODE=$?
    echo "Performing cleanup..."
    
    # Deactivate conda environment if it was activated
    if [[ -n "$CONDA_DEFAULT_ENV" ]] && [[ "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "Deactivating conda environment: $CONDA_DEFAULT_ENV"
        conda deactivate
    fi
    
    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "Job finished successfully."
    else
        echo "Job failed with exit code $EXIT_CODE."
    fi
    echo "Cleanup complete."
}
trap cleanup EXIT

# Activate bioinformatics conda environment
ENV_NAME="atac_analysis"
echo "Activating conda environment: $ENV_NAME"

# Initialize conda for bash (required for non-interactive shells)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh

if conda env list | grep -q "^${ENV_NAME} "; then
    conda activate $ENV_NAME
    echo "✓ Activated environment: $ENV_NAME"
    echo "Python version: $(python --version)"
else
    echo "ERROR: Conda environment '$ENV_NAME' not found!"
    echo "Please run: bash setup_bioinformatics_env.sh"
    exit 1
fi

# --- Configuration ---
# Tools - using conda environment tools
TRIMMOMATIC="/beegfs/scratch/ric.sessa/kubacki.michal/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
SAMTOOLS="samtools"
MACS2="macs2"
BWA="bwa"
BOWTIE2="bowtie2"
PICARD="/beegfs/scratch/ric.sessa/kubacki.michal/tools/picard.jar"

# Verify tools are available in conda environment
echo "Verifying tools in conda environment..."
if command -v bwa &> /dev/null; then
    echo "✓ BWA available: $(which bwa)"
else
    echo "✗ BWA not found in environment"
    exit 1
fi

if command -v samtools &> /dev/null; then
    echo "✓ samtools available: $(which samtools)"
else
    echo "✗ samtools not found in environment"
    exit 1
fi

if command -v macs2 &> /dev/null; then
    echo "✓ MACS2 available: $(which macs2)"
else
    echo "! MACS2 not found (peak calling will be skipped)"
fi

# Reference Genome
REFERENCE_FASTA="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-cellranger-arc-GRCm39-2024-A/fasta/genome.fa"
GENOME_SIZE="mm" # Mouse genome size for MACS2
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/GRCm39.chrom.sizes"

# Input/Output Directories - Updated for EMX1
SOURCE_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/emx1"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/bulk_atac_emx1_output"
TMP_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/tmp"

# Sample configuration - EMX1 sample
SAMPLES=("R26-Emx1-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Processing bulk ATAC-seq for EMX1 sample: $SAMPLE"
echo "Date: $(date)"
echo "=========================================="

# Create directories
mkdir -p "$OUTPUT_DIR/$SAMPLE" "$TMP_DIR" "logs"

# Step 1: Prepare genomic reads (R1 and R3)
echo "Step 1: Preparing genomic reads..."
R1_FILE="$SOURCE_DATA_DIR/${SAMPLE}_R1_001.fastq.gz"
R3_FILE="$SOURCE_DATA_DIR/${SAMPLE}_R3_001.fastq.gz"

# Verify input files exist
for file in "$R1_FILE" "$R3_FILE"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

echo "Input files verified:"
echo "  R1 (genomic): $R1_FILE"
echo "  R3 (genomic): $R3_FILE"

# Step 2: Quality trimming with Trimmomatic
echo "Step 2: Quality trimming..."
TRIMMED_DIR="$OUTPUT_DIR/$SAMPLE/trimmed"
mkdir -p "$TRIMMED_DIR"

if [[ -f "$TRIMMOMATIC" ]]; then
    java -jar "$TRIMMOMATIC" PE -threads $SLURM_CPUS_PER_TASK \
        "$R1_FILE" "$R3_FILE" \
        "$TRIMMED_DIR/${SAMPLE}_R1_paired.fq.gz" "$TRIMMED_DIR/${SAMPLE}_R1_unpaired.fq.gz" \
        "$TRIMMED_DIR/${SAMPLE}_R3_paired.fq.gz" "$TRIMMED_DIR/${SAMPLE}_R3_unpaired.fq.gz" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
else
    echo "WARNING: Trimmomatic not found, proceeding without trimming"
    TRIMMED_DIR="$SOURCE_DATA_DIR"
    # Use original files
    ln -sf "$R1_FILE" "$TRIMMED_DIR/${SAMPLE}_R1_paired.fq.gz"
    ln -sf "$R3_FILE" "$TRIMMED_DIR/${SAMPLE}_R3_paired.fq.gz"
fi

# Step 3: Alignment
echo "Step 3: Setting up alignment..."
ALIGN_DIR="$OUTPUT_DIR/$SAMPLE/alignment"
mkdir -p "$ALIGN_DIR"

# Check reference genome
if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference genome FASTA not found at $REFERENCE_FASTA"
    exit 1
fi

# Use BWA for alignment
if [[ -n "$BWA" ]]; then
    echo "Using BWA for alignment..."
    
    # Create BWA index if it doesn't exist
    if [[ ! -f "${REFERENCE_FASTA}.bwt" ]]; then
        echo "Creating BWA index (this may take 10-20 minutes)..."
        $BWA index "$REFERENCE_FASTA"
        echo "BWA index created successfully"
    fi
    
    # Align with BWA MEM
    echo "Starting BWA alignment..."
    $BWA mem -t $SLURM_CPUS_PER_TASK -M "$REFERENCE_FASTA" \
        "$TRIMMED_DIR/${SAMPLE}_R1_paired.fq.gz" \
        "$TRIMMED_DIR/${SAMPLE}_R3_paired.fq.gz" \
        > "$ALIGN_DIR/${SAMPLE}_aligned.sam"
else
    echo "ERROR: BWA not found in conda environment"
    exit 1
fi

echo "Alignment completed successfully for $SAMPLE"

# Step 4: Convert to BAM and sort
echo "Step 4: Converting to BAM and sorting..."
$SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -bS "$ALIGN_DIR/${SAMPLE}_aligned.sam" > "$ALIGN_DIR/${SAMPLE}_aligned.bam"
$SAMTOOLS sort -@ $SLURM_CPUS_PER_TASK "$ALIGN_DIR/${SAMPLE}_aligned.bam" -o "$ALIGN_DIR/${SAMPLE}_sorted.bam"
$SAMTOOLS index "$ALIGN_DIR/${SAMPLE}_sorted.bam"

# Remove large SAM file to save space
rm "$ALIGN_DIR/${SAMPLE}_aligned.sam"

# Step 5: Filter for high-quality alignments
echo "Step 5: Filtering for high-quality alignments..."
# Filter for:
# - Properly paired reads
# - MAPQ >= 30
# - Remove duplicates
$SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -b -f 2 -q 30 "$ALIGN_DIR/${SAMPLE}_sorted.bam" > "$ALIGN_DIR/${SAMPLE}_filtered.bam"

# Mark duplicates if Picard is available
if [[ -f "$PICARD" ]]; then
    echo "Marking duplicates with Picard..."
    java -jar "$PICARD" MarkDuplicates \
        INPUT="$ALIGN_DIR/${SAMPLE}_filtered.bam" \
        OUTPUT="$ALIGN_DIR/${SAMPLE}_dedup.bam" \
        METRICS_FILE="$ALIGN_DIR/${SAMPLE}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true
    
    $SAMTOOLS index "$ALIGN_DIR/${SAMPLE}_dedup.bam"
    FINAL_BAM="$ALIGN_DIR/${SAMPLE}_dedup.bam"
else
    echo "Picard not found, skipping duplicate removal"
    FINAL_BAM="$ALIGN_DIR/${SAMPLE}_filtered.bam"
    $SAMTOOLS index "$FINAL_BAM"
fi

# Step 6: Generate alignment statistics
echo "Step 6: Generating alignment statistics..."
STATS_DIR="$OUTPUT_DIR/$SAMPLE/stats"
mkdir -p "$STATS_DIR"

$SAMTOOLS flagstat "$FINAL_BAM" > "$STATS_DIR/${SAMPLE}_flagstat.txt"
$SAMTOOLS idxstats "$FINAL_BAM" > "$STATS_DIR/${SAMPLE}_idxstats.txt"

# Step 7: Peak calling with MACS2
echo "Step 7: Calling peaks with MACS2..."
PEAKS_DIR="$OUTPUT_DIR/$SAMPLE/peaks"
mkdir -p "$PEAKS_DIR"

# Call peaks for ATAC-seq data (if MACS2 is available)
if command -v macs2 &> /dev/null; then
    echo "Running MACS2 peak calling..."
    $MACS2 callpeak \
        -t "$FINAL_BAM" \
        -f BAMPE \
        -g "$GENOME_SIZE" \
        -n "$SAMPLE" \
        --outdir "$PEAKS_DIR" \
        --shift -75 \
        --extsize 150 \
        --nomodel \
        --broad \
        --broad-cutoff 0.1 \
        --keep-dup all
    echo "Peak calling completed"
else
    echo "MACS2 not available - skipping peak calling"
    echo "You can run peak calling manually later:"
    echo "  conda activate atac_analysis"
    echo "  macs2 callpeak -t $FINAL_BAM -f BAMPE -g mm -n $SAMPLE --outdir $PEAKS_DIR"
fi

# Step 8: Generate summary report
echo "Step 8: Generating summary report..."
REPORT_FILE="$OUTPUT_DIR/$SAMPLE/${SAMPLE}_processing_report.txt"

cat > "$REPORT_FILE" << EOF
Bulk ATAC-seq Processing Report - EMX1 Sample
==============================================
Sample: $SAMPLE
Processing Date: $(date)
Sample Type: EMX1 Mutant Adult

Input Files:
- R1 (genomic): $R1_FILE
- R3 (genomic): $R3_FILE

Output Files:
- Final BAM: $FINAL_BAM
- Peaks (narrowPeak): $PEAKS_DIR/${SAMPLE}_peaks.narrowPeak
- Peaks (broadPeak): $PEAKS_DIR/${SAMPLE}_peaks.broadPeak
- Statistics: $STATS_DIR/

Alignment Statistics:
$(cat $STATS_DIR/${SAMPLE}_flagstat.txt)

Peak Calling Results:
Total peaks called: $(wc -l < $PEAKS_DIR/${SAMPLE}_peaks.narrowPeak || echo "N/A")
Broad peaks called: $(wc -l < $PEAKS_DIR/${SAMPLE}_peaks.broadPeak || echo "N/A")

Processing completed successfully!
EOF

echo "=========================================="
echo "Bulk ATAC-seq processing completed for EMX1 sample: $SAMPLE"
echo "Results available at: $OUTPUT_DIR/$SAMPLE"
echo "Summary report: $REPORT_FILE"
echo "=========================================="

# Optional: Generate browser tracks
if [[ -f "$CHROM_SIZES" ]]; then
    echo "Generating browser tracks..."
    TRACKS_DIR="$OUTPUT_DIR/$SAMPLE/tracks"
    mkdir -p "$TRACKS_DIR"
    
    # Convert BAM to bedGraph
    $SAMTOOLS depth "$FINAL_BAM" | \
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > "$TRACKS_DIR/${SAMPLE}_coverage.bedgraph"
    
    echo "Browser tracks generated at: $TRACKS_DIR"
fi

echo "Processing pipeline completed for EMX1 sample: $SAMPLE"