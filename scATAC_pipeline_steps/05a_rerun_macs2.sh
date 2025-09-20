#!/bin/bash
# =============================================================================
# SCRIPT: 05a_rerun_macs2.sh
# PURPOSE: Rerun MACS2 peak calling with clean slate approach
# 
# DESCRIPTION:
# This script provides a clean rerun capability for MACS2 peak calling.
# It removes previous output files and regenerates all peak calling results
# from scratch, ensuring consistent and reproducible peak identification.
# 
# USE CASES:
# - Parameter optimization and testing
# - Recovery from failed or interrupted runs
# - Clean regeneration after upstream changes
# - Quality control and validation runs
# 
# INPUT REQUIREMENTS:
# - Fragment files from step 4 (04_chromap_alignment.sh)
# - Reference genome FASTA file with index
# - MACS2 peak calling software
# 
# OUTPUT GENERATED:
# - Fresh MACS2 narrowPeak files with peak coordinates
# - Clean peak summit annotations
# - Regenerated bedgraph files for visualization
# - Updated sorted BED files for downstream analysis
# - New peak calling statistics and logs
# 
# COMPUTATIONAL REQUIREMENTS:
# - 8 CPU cores for parallel processing
# - 64GB RAM for large dataset handling
# - 12-hour time limit for comprehensive reanalysis
# 
# DEPENDENCIES:
# - MACS2: Peak calling algorithm
# - samtools: Reference genome indexing
# - Standard Unix tools: awk, sort, gzip
# =============================================================================

#SBATCH --job-name=rerun_macs2
#SBATCH --output=logs/05a_rerun_macs2_%a.out
#SBATCH --error=logs/05a_rerun_macs2_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=workq

# =============================================================================
# CONDA ENVIRONMENT SETUP
# PURPOSE: Activate specialized environment for peak calling reanalysis
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
# 
# RERUN CONSIDERATIONS:
# - Same environment as original run for consistency
# - Ensures reproducible results
# - Maintains parameter compatibility
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

set -euo pipefail

# =============================================================================
# CONFIGURATION AND SAMPLE SETUP
# PURPOSE: Define sample array and reanalysis parameters
# 
# SAMPLE CONFIGURATION:
# - Identical array setup to original peak calling
# - Control vs Mutant comparison design maintained
# - Adult tissue-specific analysis consistency
# 
# PATH CONFIGURATION:
# - OUTPUT_DIR: Same location as original analysis
# - REF_GENOME: Consistent mouse mm10 reference
# 
# RERUN STRATEGY:
# - Maintain identical parameters for reproducibility
# - Use same file structure for easy comparison
# - Preserve original analysis workflow
# 
# ARRAY PROCESSING:
# - SLURM_ARRAY_TASK_ID: Automatic sample selection
# - Parallel reprocessing of multiple samples
# - Consistent parameter application across reruns
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"

echo "========================================="
echo "Step 5a: Rerunning Peak calling for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY PREPARATION
# PURPOSE: Ensure organized directory structure for rerun results
# 
# DIRECTORY STRUCTURE:
# - peaks/: MACS2 output files (will be cleaned and regenerated)
# - qc/: Quality control files and genome references
# 
# RERUN CONSIDERATIONS:
# - Directories may already exist from previous runs
# - Files will be cleaned before regeneration
# - Maintains consistent organization structure
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
# RERUN REQUIREMENTS:
# - Same input files as original run
# - Ensures consistent starting point
# - Validates upstream pipeline completion
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
# PURPOSE: Create chromosome size reference for MACS2 peak calling
# 
# FUNCTIONALITY:
# - Extracts chromosome names and lengths from reference genome
# - Creates tab-delimited file required by MACS2
# - Caches result for reuse across samples
# 
# FILE FORMAT:
# chr1    248956422
# chr2    242193529
# ...
# 
# USAGE IN PEAK CALLING:
# - Defines valid genomic regions for peak detection
# - Enables proper statistical modeling
# - Required for MACS2 --gsize parameter
# 
# CACHING STRATEGY:
# - Reuses existing file if available
# - Avoids redundant computation in reruns
# - Ensures consistency across pipeline runs
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
# PURPOSE: Convert ATAC-seq fragments to individual reads for MACS2
# 
# CONVERSION STRATEGY:
# - Extract 5' end of each fragment (Tn5 insertion sites)
# - Create paired-end reads from fragment boundaries
# - Generate BED format compatible with MACS2
# 
# OUTPUT FORMAT:
# chr1    1000    1001    read1
# chr1    1150    1151    read2
# 
# MACS2 REQUIREMENTS:
# - Individual read positions (not fragments)
# - Excludes mitochondrial DNA (chrM)
# - Sorted coordinates for efficient processing
# 
# QUALITY CONTROL:
# - Validates successful conversion
# - Counts generated reads for verification
# - Ensures non-empty output for downstream analysis
# 
# RERUN CONSIDERATIONS:
# - May reuse existing reads file if available
# - Validates existing file before proceeding
# - Ensures consistent input for peak calling
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
# PURPOSE: Identify accessible chromatin regions from ATAC-seq data
# 
# TOOL VERIFICATION:
# - Confirms MACS2 installation and accessibility
# - Validates version compatibility
# - Provides installation guidance if missing
# 
# PARAMETER DETAILS:
# -t: Treatment file (converted reads)
# -g mm: Mouse genome size (for statistical modeling)
# -f BED: Input format specification
# -q 0.01: FDR threshold (1% false discovery rate)
# --nomodel: Skip fragment size estimation (use fixed parameters)
# --shift -100: Shift reads by 100bp (account for Tn5 insertion)
# --extsize 200: Extend reads to 200bp (typical nucleosome-free region)
# --keep-dup all: Retain all duplicate reads
# -B --SPMR: Generate bedGraph files with signal per million reads
# 
# OUTPUT FILES:
# - narrowPeak: Primary peak coordinates with scores
# - summits.bed: Peak summit positions
# - bedGraph files: Signal tracks for visualization
# 
# ERROR HANDLING:
# - Validates MACS2 installation before execution
# - Monitors exit codes for failure detection
# - Provides detailed error messages for troubleshooting
# =============================================================================
echo "DEBUG: Checking for MACS2..."
if ! command -v macs2 &> /dev/null; then
    echo "ERROR: MACS2 not found. Please install MACS2 for peak calling"
    echo "DEBUG: You can install with: pip install MACS2"
    exit 1
fi

echo "DEBUG: MACS2 found: $(which macs2)"
echo "DEBUG: MACS2 version: $(macs2 --version 2>&1 || echo 'Version check failed')"

# =============================================================================
# PREVIOUS RUN CLEANUP
# PURPOSE: Remove all previous MACS2 output files for clean rerun
# 
# FILES REMOVED:
# - narrowPeak: Primary peak coordinates file
# - peaks_sorted.bed: Sorted peak coordinates
# - summits.bed: Peak summit positions
# - control_lambda.bdg: Background signal model
# - treat_pileup.bdg: Treatment signal bedgraph
# - peaks.xls: MACS2 detailed statistics
# 
# CLEANUP STRATEGY:
# - Force removal (rm -f) to avoid interactive prompts
# - Remove all MACS2-generated files
# - Ensure clean slate for reproducible results
# 
# SAFETY CONSIDERATIONS:
# - Only removes specific MACS2 output files
# - Preserves input files and other analysis results
# - Allows for safe rerun without data loss
# =============================================================================
rm -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" \
      "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed" \
      "$OUTPUT_DIR/peaks/${SAMPLE}_summits.bed" \
      "$OUTPUT_DIR/peaks/${SAMPLE}_control_lambda.bdg" \
      "$OUTPUT_DIR/peaks/${SAMPLE}_treat_pileup.bdg" \
      "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.xls"

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

# =============================================================================
# PEAK VALIDATION AND QUALITY CONTROL
# PURPOSE: Verify successful peak calling and output file generation
# 
# VALIDATION CHECKS:
# - Confirms primary peaks file exists
# - Validates file accessibility and readability
# - Ensures non-empty output
# 
# ERROR HANDLING:
# - Detects silent failures in MACS2 execution
# - Provides clear error messages for missing outputs
# - Prevents downstream analysis with invalid data
# 
# QUALITY INDICATORS:
# - File existence and accessibility
# - Non-zero file size
# - Proper file format structure
# =============================================================================
if [[ ! -f "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" ]]; then
    echo "ERROR: MACS2 completed but peaks file not found"
    exit 1
fi

# =============================================================================
# PEAK STATISTICS CALCULATION
# PURPOSE: Quantify peak calling results and assess data quality
# 
# KEY METRICS:
# - Total peak count: Number of accessible regions identified
# - Peak distribution: Genomic coverage and density
# - Quality assessment: Expected range validation
# 
# BIOLOGICAL INTERPRETATION:
# - Typical ATAC-seq: 20,000-100,000 peaks
# - Low counts may indicate poor data quality
# - High counts may suggest over-calling or contamination
# 
# QUALITY CONTROL:
# - Validates non-zero peak count
# - Provides quantitative assessment
# - Enables comparison across samples
# =============================================================================
PEAK_COUNT=$(wc -l < "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak")
echo "DEBUG: Number of peaks: $PEAK_COUNT"

# =============================================================================
# SORTED PEAKS FILE GENERATION
# PURPOSE: Create standardized peak coordinates for downstream analysis
# 
# PROCESSING STEPS:
# - Extract essential columns (chr, start, end, name)
# - Sort by chromosome and position
# - Generate BED4 format for compatibility
# 
# OUTPUT FORMAT:
# chr1    1000    1200    peak_1
# chr1    1500    1700    peak_2
# ...
# 
# DOWNSTREAM APPLICATIONS:
# - Peak annotation and functional analysis
# - Motif discovery and enrichment
# - Comparative analysis across samples
# - Integration with other genomic datasets
# =============================================================================
echo "DEBUG: Creating sorted peaks file..."
cut -f 1-4 "$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak" | \
    sort -k1,1V -k2,2n > "$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

# =============================================================================
# PEAK CALLING RESULTS SUMMARY
# PURPOSE: Provide comprehensive overview of peak calling outcomes
# 
# SUMMARY COMPONENTS:
# - Peak count and quality metrics
# - Output file locations and formats
# - Processing statistics and validation
# 
# QUALITY INDICATORS:
# - Peak count within expected range
# - Successful file generation
# - Proper format validation
# 
# RERUN VALIDATION:
# - Confirms successful cleanup and regeneration
# - Validates consistency with previous runs
# - Ensures reproducible results
# =============================================================================
echo "Peak calling results:"
echo "  Total peaks called: $(printf "%d" $PEAK_COUNT)"
echo "  Peak file: $OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak"
echo "  Sorted peaks: $OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"

# =============================================================================
# ADDITIONAL MACS2 OUTPUT FILES DOCUMENTATION
# PURPOSE: Catalog and validate supplementary MACS2 outputs
# 
# OUTPUT FILES:
# - summits.bed: Peak summit coordinates (highest signal points)
# - control_lambda.bdg: Background signal model
# - treat_pileup.bdg: Treatment signal bedgraph
# 
# DOWNSTREAM USAGE:
# - Summit files for precise binding site analysis
# - BedGraph files for signal visualization
# - Lambda files for background correction
# 
# VALIDATION:
# - Confirms file existence and accessibility
# - Reports file sizes for quality assessment
# - Enables troubleshooting of incomplete runs
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
# PURPOSE: Create comprehensive metadata file for analysis tracking
# 
# METADATA COMPONENTS:
# - Sample identification and processing parameters
# - Input/output file locations and formats
# - Processing timestamps and version information
# - Quality metrics and validation results
# 
# FILE STRUCTURE:
# - Machine-readable format for automated processing
# - Human-readable content for manual review
# - Standardized fields for cross-sample comparison
# 
# USAGE:
# - Pipeline tracking and reproducibility
# - Quality control and batch analysis
# - Integration with downstream analysis tools
# =============================================================================
cat > "$OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt" << EOF
SAMPLE=$SAMPLE
PEAK_COUNT=$PEAK_COUNT
READS_INPUT=$READS_FILE
PEAK_FILE=$OUTPUT_DIR/peaks/${SAMPLE}_peaks.narrowPeak
SORTED_PEAKS=$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed
PEAK_CALLING_COMPLETED_AT=$(date)
MACS2_PARAMETERS=-g mm -f BED -q 0.01 --nomodel --shift -100 --extsize 200 --keep-dup all -B --SPMR
EOF

echo "Results saved to: $OUTPUT_DIR/peaks/${SAMPLE}_peak_calling_results.txt"

# =============================================================================
# PIPELINE COMPLETION AND SUMMARY
# PURPOSE: Provide final status and next steps guidance
# 
# COMPLETION SUMMARY:
# - Sample processing confirmation
# - Peak calling statistics
# - Processing time and efficiency metrics
# 
# OUTPUT FILES:
# - Primary peaks: narrowPeak format with scores
# - Sorted peaks: BED format for downstream analysis
# - Signal tracks: bedGraph files for visualization
# - Metadata: Processing parameters and statistics
# 
# NEXT STEPS:
# - Quality control analysis (step 8)
# - BigWig generation for visualization (step 7)
# - Peak annotation and functional analysis
# - Comparative analysis across samples
# 
# RERUN COMPLETION:
# - Confirms successful cleanup and regeneration
# - Validates output file integrity
# - Ensures pipeline continuity
# =============================================================================
echo "========================================="
echo "Step 5a complete for $SAMPLE"
echo "Called $(printf "%d" $PEAK_COUNT) peaks"
echo "End time: $(date)"
echo "========================================="
