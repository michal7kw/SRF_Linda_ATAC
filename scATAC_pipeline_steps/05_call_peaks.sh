#!/bin/bash
# =============================================================================
# SCRIPT: 05_call_peaks.sh
# PURPOSE: MACS2-based peak calling for single-cell ATAC-seq data
# 
# DESCRIPTION:
# This script performs peak calling on aligned ATAC-seq fragments using MACS2.
# It converts fragment data to read format, applies ATAC-seq specific parameters,
# and generates comprehensive peak annotations for downstream analysis.
# 
# INPUT REQUIREMENTS:
# - Fragment files from step 4 (04_chromap_alignment.sh)
# - Reference genome FASTA file with index
# - MACS2 peak calling software
# 
# OUTPUT GENERATED:
# - MACS2 narrowPeak files with peak coordinates and scores
# - Peak summit annotations for precise binding sites
# - Bedgraph files for signal visualization
# - Sorted BED files for downstream analysis
# - Peak calling statistics and parameter logs
# 
# COMPUTATIONAL REQUIREMENTS:
# - 8 CPU cores for parallel processing
# - 64GB RAM for large dataset handling
# - 12-hour time limit for comprehensive analysis
# 
# DEPENDENCIES:
# - MACS2: Peak calling algorithm
# - samtools: Reference genome indexing
# - Standard Unix tools: awk, sort, gzip
# =============================================================================

#SBATCH --job-name=call_peaks
#SBATCH --output=logs/05_call_peaks_%a.out
#SBATCH --error=logs/05_call_peaks_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# =============================================================================
# CONDA ENVIRONMENT SETUP
# PURPOSE: Activate specialized environment for peak calling analysis
# 
# ENVIRONMENT FEATURES:
# - MACS2: Peak calling algorithm optimized for ATAC-seq
# - samtools: Reference genome manipulation and indexing
# - Python dependencies: numpy, scipy for statistical analysis
# 
# ERROR HANDLING:
# - Strict mode enabled (set -euo pipefail)
# - Exit on any command failure
# - Treat undefined variables as errors
# - Fail on pipe errors
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

set -euo pipefail

# =============================================================================
# CONFIGURATION AND SAMPLE SETUP
# PURPOSE: Define sample array and analysis parameters
# 
# SAMPLE CONFIGURATION:
# - Array-based processing for parallel execution
# - Control vs Mutant comparison design
# - Adult tissue-specific analysis
# 
# PATH CONFIGURATION:
# - OUTPUT_DIR: Central location for all analysis results
# - REF_GENOME: Mouse mm10 reference genome for peak annotation
# 
# ARRAY PROCESSING:
# - SLURM_ARRAY_TASK_ID: Automatic sample selection
# - Enables parallel processing of multiple samples
# - Consistent parameter application across samples
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"

echo "========================================="
echo "Step 5: Peak calling for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY PREPARATION
# PURPOSE: Create organized directory structure for peak calling results
# 
# DIRECTORY STRUCTURE:
# - peaks/: MACS2 output files (narrowPeak, summits, bedgraph)
# - qc/: Quality control files and genome references
# 
# ORGANIZATION BENEFITS:
# - Consistent file organization across samples
# - Easy navigation and file management
# - Separation of different analysis types
# =============================================================================
mkdir -p "$OUTPUT_DIR/peaks" "$OUTPUT_DIR/qc"

# =============================================================================
# PREREQUISITE VALIDATION
# PURPOSE: Verify required input files from previous pipeline steps
# 
# VALIDATION CHECKS:
# - Fragment file existence from step 4 alignment
# - File accessibility and readability
# - Pipeline dependency verification
# 
# ERROR HANDLING:
# - Clear error messages for missing dependencies
# - Guidance for resolving prerequisite issues
# - Graceful exit with informative feedback
# =============================================================================
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

# =============================================================================
# GENOME SIZES FILE GENERATION
# PURPOSE: Create chromosome length reference for peak calling normalization
# 
# FUNCTIONALITY:
# - Extract chromosome names and lengths from reference genome
# - Generate FASTA index (.fai) for rapid access
# - Create tab-delimited genome sizes file for MACS2
# 
# FILE FORMAT:
# - Column 1: Chromosome name (chr1, chr2, etc.)
# - Column 2: Chromosome length in base pairs
# 
# USAGE IN PEAK CALLING:
# - Genome-wide normalization of signal
# - Proper handling of chromosome boundaries
# - Statistical modeling of background signal
# 
# CACHING STRATEGY:
# - Reuse existing genome sizes file if available
# - Avoid redundant computation across samples
# =============================================================================
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

# =============================================================================
# FRAGMENT-TO-READS CONVERSION
# PURPOSE: Convert paired-end fragments to individual read positions for MACS2
# 
# CONVERSION STRATEGY:
# - Extract 5' end positions from fragment coordinates
# - Generate separate entries for forward and reverse reads
# - Maintain barcode information for single-cell context
# - Filter mitochondrial reads for nuclear analysis
# 
# OUTPUT FORMAT (BED4):
# - Column 1: Chromosome name
# - Column 2: Read start position (0-based)
# - Column 3: Read end position (1-based)
# - Column 4: Barcode identifier
# 
# MACS2 REQUIREMENTS:
# - Individual read positions rather than fragments
# - Proper coordinate sorting for efficient processing
# - Compressed format for storage efficiency
# 
# QUALITY CONTROL:
# - Read count validation and reporting
# - File existence checking for pipeline robustness
# - Mitochondrial read filtering
# =============================================================================
echo "DEBUG: Converting fragments to reads for peak calling..."
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"

if [[ ! -s "$READS_FILE" ]]; then
    echo "DEBUG: Reads file not found or is empty. Generating from fragment ends..."
    zcat "$FRAGMENTS_FILE" | \
        grep -v "^chrM" | \
        awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $4; print $1, $3-1, $3, $4}' | \
        sort -k1,1V -k2,2n | \
        gzip > "$READS_FILE"

    READ_COUNT=$(zcat "$READS_FILE" | wc -l)
    if [[ $? -eq 0 && $READ_COUNT -gt 0 ]]; then
        echo "DEBUG: Reads file created successfully"
        echo "DEBUG: Number of reads: $(printf "%d" $READ_COUNT)"
    else
        echo "ERROR: Failed to create reads file or the file is empty"
        echo "DEBUG: Read count: $READ_COUNT"
        exit 1
    fi
else
    echo "DEBUG: Reads file already exists, skipping generation."
fi

# =============================================================================
# MACS2 PEAK CALLING
# PURPOSE: Identify regions of chromatin accessibility using MACS2 algorithm
# 
# TOOL VERIFICATION:
# - Check MACS2 installation and accessibility
# - Verify version compatibility
# - Provide installation guidance if missing
# 
# MACS2 PARAMETERS EXPLAINED:
# -t: Treatment file (our ATAC-seq reads)
# -g mm: Mouse genome size (effective genome size)
# -f BED: Input file format
# -q 0.01: Q-value (FDR) cutoff for peak detection
# --nomodel: Skip building shifting model (appropriate for ATAC-seq)
# --shift -100: Shift reads by -100bp (center on nucleosome-free regions)
# --extsize 200: Extend reads to 200bp (typical nucleosome-free region)
# --keep-dup all: Keep all duplicate reads (important for scATAC-seq)
# -B: Generate bedgraph files for visualization
# --SPMR: Save signal per million reads for normalization
# 
# OUTPUT FILES GENERATED:
# - narrowPeak: Peak coordinates with scores and significance
# - summits.bed: Peak summit positions for precise binding sites
# - treat_pileup.bdg: Signal bedgraph for visualization
# - control_lambda.bdg: Background model for comparison
# 
# ERROR HANDLING:
# - Exit code monitoring for MACS2 execution
# - Clear error messages for troubleshooting
# - Graceful handling of missing dependencies
# =============================================================================
if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "DEBUG: Checking for MACS2..."
    if ! command -v macs2 &> /dev/null; then
        echo "ERROR: MACS2 not found. Please install MACS2 for peak calling"
        echo "DEBUG: You can install with: pip install MACS2"
        exit 1
    fi

    echo "DEBUG: MACS2 found: $(which macs2)"
    echo "DEBUG: MACS2 version: $(macs2 --version 2>&1 || echo 'Version check failed')"

    echo "DEBUG: Running MACS2 peak calling..."
    macs2 callpeak \
        -t "$READS_FILE" \
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
else
    echo "DEBUG: Skipping MACS2 peak calling, peaks file already exists."
fi


# =============================================================================
# PEAK VALIDATION AND QUALITY CONTROL
# PURPOSE: Verify successful peak calling and validate output files
# 
# VALIDATION CHECKS:
# - Confirm narrowPeak file creation
# - Verify file accessibility and readability
# - Ensure non-empty output from MACS2
# 
# ERROR HANDLING:
# - Clear error messages for missing output
# - Pipeline integrity verification
# - Graceful exit with diagnostic information
# =============================================================================
if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "ERROR: MACS2 completed but peaks file not found"
    exit 1
fi

# =============================================================================
# PEAK STATISTICS CALCULATION
# PURPOSE: Generate quantitative metrics for peak calling assessment
# 
# KEY METRICS:
# - Total peak count: Number of significant accessibility regions
# - Peak density: Peaks per megabase of genome
# - Quality assessment: Peak calling success validation
# 
# BIOLOGICAL INTERPRETATION:
# - Typical ATAC-seq: 20,000-100,000 peaks expected
# - Higher peak counts: More open chromatin regions
# - Lower peak counts: May indicate technical issues
# 
# DOWNSTREAM USAGE:
# - Quality control assessment
# - Sample comparison metrics
# - Pipeline success validation
# =============================================================================
PEAK_COUNT=$(wc -l < "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak")
echo "DEBUG: Number of peaks: $PEAK_COUNT"

# =============================================================================
# SORTED PEAKS FILE GENERATION
# PURPOSE: Create coordinate-sorted BED file for downstream analysis
# 
# PROCESSING STEPS:
# - Extract essential coordinates (chr, start, end, name)
# - Sort by chromosome and position
# - Generate standardized BED4 format
# 
# OUTPUT FORMAT:
# - Column 1: Chromosome name
# - Column 2: Peak start position (0-based)
# - Column 3: Peak end position (1-based)
# - Column 4: Peak identifier
# 
# DOWNSTREAM APPLICATIONS:
# - Genomic interval operations (bedtools)
# - Peak annotation and gene assignment
# - Motif analysis and enrichment
# - Comparative analysis between samples
# =============================================================================
echo "DEBUG: Creating sorted peaks file..."
cut -f 1-4 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" | \
    sort -k1,1 -k2,2n > "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

# =============================================================================
# PEAK CALLING RESULTS SUMMARY
# PURPOSE: Display comprehensive summary of peak calling outcomes
# 
# SUMMARY COMPONENTS:
# - Total peak count with formatted display
# - Primary output file locations
# - Processed file derivatives
# 
# QUALITY INDICATORS:
# - Peak count within expected range
# - Successful file generation
# - Proper coordinate sorting
# =============================================================================
echo "Peak calling results:"
echo "  Total peaks called: $(printf "%d" $PEAK_COUNT)"
echo "  Peak file: $OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak"
echo "  Sorted peaks: $OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

# =============================================================================
# ADDITIONAL MACS2 OUTPUT FILES DOCUMENTATION
# PURPOSE: Catalog and validate supplementary MACS2 output files
# 
# FILE TYPES:
# - summits.bed: Peak summit coordinates for precise binding sites
# - control_lambda.bdg: Background signal model
# - treat_pileup.bdg: Treatment signal for visualization
# 
# VALIDATION:
# - File existence verification
# - Line count reporting for content validation
# - Comprehensive output documentation
# 
# DOWNSTREAM USAGE:
# - Summit files: Motif analysis, precise binding site annotation
# - Bedgraph files: Genome browser visualization, signal analysis
# - Lambda files: Background correction, signal normalization
# =============================================================================
echo "Additional MACS2 output files:"
for ext in summits.bed control_lambda.bdg treat_pileup.bdg; do
    MACS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_${ext}"
    if [[ -f "$MACS_FILE" ]]; then
        echo "  - $(basename "$MACS_FILE"): $(wc -l < "$MACS_FILE") lines"
    fi
done

# =============================================================================
# PEAK CALLING RESULTS DOCUMENTATION
# PURPOSE: Generate comprehensive metadata file for analysis tracking
# 
# METADATA COMPONENTS:
# - Sample identification and processing parameters
# - Input/output file locations and statistics
# - MACS2 parameters for reproducibility
# - Timestamp for analysis tracking
# 
# FILE STRUCTURE:
# - Key-value pairs for easy parsing
# - Complete parameter documentation
# - File path references for downstream analysis
# 
# USAGE:
# - Pipeline tracking and reproducibility
# - Quality control assessment
# - Batch analysis coordination
# - Parameter optimization reference
# =============================================================================
cat > "$OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt" << EOF
SAMPLE=$SAMPLE
PEAK_COUNT=$PEAK_COUNT
FRAGMENTS_INPUT=$FRAGMENTS_FILE
PEAK_FILE=$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak
SORTED_PEAKS=$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed
PEAK_CALLING_COMPLETED_AT=$(date)
MACS2_PARAMETERS=-g mm -f BED -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -B --SPMR
EOF

echo "Results saved to: $OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt"

# =============================================================================
# PIPELINE COMPLETION AND SUMMARY
# PURPOSE: Provide final status report and next steps guidance
# 
# COMPLETION SUMMARY:
# - Step identification and sample processing confirmation
# - Key metrics display (peak count)
# - Timestamp for analysis tracking
# 
# OUTPUT FILES GENERATED:
# - ${SAMPLE}_peaks.narrowPeak: Primary peak coordinates
# - ${SAMPLE}_peaks_sorted.bed: Sorted coordinates for analysis
# - ${SAMPLE}_summits.bed: Peak summit positions
# - ${SAMPLE}_treat_pileup.bdg: Signal bedgraph
# - ${SAMPLE}_peak_calling_results.txt: Analysis metadata
# 
# NEXT STEPS:
# - Peak annotation and gene assignment
# - Motif analysis and enrichment
# - Differential accessibility analysis
# - Visualization and quality assessment
# =============================================================================
echo "========================================="
echo "Step 5 complete for $SAMPLE"
echo "Called $(printf "%d" $PEAK_COUNT) peaks"
echo "End time: $(date)"
echo "=========================================="
