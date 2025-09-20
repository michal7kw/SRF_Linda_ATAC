#!/bin/bash
# =============================================================================
# SCRIPT: 04b_alignment_qc.sh
# PURPOSE: Comprehensive quality control analysis for chromap alignment results
# 
# DESCRIPTION:
# This script performs detailed quality control analysis on chromap alignment
# outputs, generating comprehensive metrics and visualizations to assess
# alignment quality, fragment characteristics, barcode performance, and
# overall data suitability for downstream peak calling analysis.
# 
# KEY OPERATIONS:
# 1. Parse chromap alignment statistics and mapping rates
# 2. Analyze fragment distribution and length characteristics
# 3. Assess chromosome distribution and mitochondrial content
# 4. Evaluate barcode quality and cell detection metrics
# 5. Generate comprehensive QC reports and recommendations
# 6. Create visualization plots for quality assessment
# 7. Provide data quality recommendations for peak calling
# 
# INPUT REQUIREMENTS:
# - Fragment files from chromap alignment (step 4)
# - Chromap alignment logs and statistics
# - Read files in BED format
# 
# OUTPUT:
# - Comprehensive QC summary reports
# - Fragment and barcode quality statistics
# - Quality assessment plots (if R available)
# - Data quality recommendations for next steps
# 
# COMPUTATIONAL REQUIREMENTS:
# - Moderate memory (64GB) for large fragment file processing
# - Multi-threading (8 cores) for efficient analysis
# - Shorter runtime (4h) focused on analysis rather than alignment
# =============================================================================
#SBATCH --job-name=alignment_qc
#SBATCH --output=logs/04b_alignment_qc_%a.out
#SBATCH --error=logs/04b_alignment_qc_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=workq

# =============================================================================
# ENVIRONMENT SETUP
# PURPOSE: Initialize conda environment and error handling for QC analysis
# 
# OPERATIONS:
# - Activates alignment_two conda environment
# - Enables strict error handling for robust analysis
# 
# REQUIRED TOOLS:
# - bedtools: Genomic interval operations and analysis
# - samtools: SAM/BAM file processing utilities
# - R/Rscript: Statistical analysis and plotting (optional)
# - Standard Unix utilities: awk, sort, uniq, bc for data processing
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# QC ANALYSIS CONFIGURATION
# PURPOSE: Define samples and output directories for quality control analysis
# 
# SAMPLE PROCESSING:
# - Uses SLURM array for parallel QC analysis of multiple samples
# - Each job analyzes one sample's alignment results independently
# 
# DIRECTORY STRUCTURE:
# - OUTPUT_DIR: Base directory containing alignment results
# - alignment_qc/: Dedicated subdirectory for QC outputs and reports
# 
# QC SCOPE:
# - Post-alignment analysis focusing on data quality assessment
# - Preparation for informed decision-making on peak calling parameters
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 4b: Alignment Quality Control for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# QC OUTPUT DIRECTORY PREPARATION
# PURPOSE: Create dedicated directory structure for QC analysis outputs
# 
# DIRECTORY ORGANIZATION:
# - alignment_qc/: Main QC analysis directory
# - Separate files for different QC metrics and analyses
# - Organized structure for easy result interpretation and sharing
# =============================================================================
mkdir -p "$OUTPUT_DIR/alignment_qc"

# =============================================================================
# PREREQUISITE VALIDATION
# PURPOSE: Verify that all required alignment output files exist for QC analysis
# 
# REQUIRED INPUT FILES:
# - FRAGMENTS_FILE: Compressed fragment file from chromap alignment
# - READS_FILE: BED format read file for duplicate analysis
# - CHROMAP_LOG: Alignment log file containing mapping statistics
# 
# VALIDATION STRATEGY:
# - Check file existence before proceeding with analysis
# - Provide clear error messages with guidance on missing dependencies
# - Ensure QC analysis has complete input data for accurate assessment
# 
# ERROR HANDLING:
# - Graceful failure with informative messages
# - Clear indication of which prerequisite step needs completion
# =============================================================================
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
CHROMAP_LOG="$OUTPUT_DIR/logs/${SAMPLE}_chromap.log"

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run step 4 (04_chromap_alignment.sh) first"
    exit 1
fi

echo "DEBUG: Prerequisites verified"

# =============================================================================
# TOOL AVAILABILITY VERIFICATION
# PURPOSE: Ensure all required bioinformatics tools are available for QC analysis
# 
# REQUIRED TOOLS:
# - bedtools: Essential for genomic interval operations, fragment analysis
# - samtools: Required for BAM/SAM file processing and statistics
# 
# VERIFICATION STRATEGY:
# - Check command availability before proceeding with analysis
# - Provide clear error messages for missing dependencies
# - Prevent runtime failures due to missing tools
# 
# TOOL USAGE IN QC:
# - bedtools: Fragment overlap analysis, chromosome distribution
# - samtools: Read statistics, duplicate analysis, format conversions
# =============================================================================
# Check for required tools
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found"
    exit 1
fi

echo "DEBUG: Required tools verified"

# =============================================================================
# CHROMAP ALIGNMENT STATISTICS PARSING
# PURPOSE: Extract and calculate key alignment metrics from chromap log files
# 
# EXTRACTED METRICS:
# - Total reads: Input read count for alignment
# - Mapped reads: Successfully aligned read count
# - Mapping rate: Percentage of reads successfully aligned
# - Duplicate rate: Percentage of PCR/optical duplicates detected
# 
# PARSING STRATEGY:
# - Use grep and awk to extract specific metrics from log files
# - Calculate derived metrics (mapping rate) using bc for precision
# - Handle missing log files gracefully with default values
# 
# QUALITY INDICATORS:
# - High mapping rate (>80%) indicates good alignment quality
# - Low duplicate rate (<20%) suggests good library complexity
# - These metrics inform downstream analysis parameter selection
# =============================================================================
# 1. Parse chromap alignment statistics
echo "DEBUG: Parsing chromap alignment statistics..."
if [[ -f "$CHROMAP_LOG" ]]; then
    # Extract key alignment metrics from chromap log
    cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" << EOF
# Chromap Alignment Statistics for $SAMPLE
# Generated: $(date)

EOF
    
    # Extract total reads
    TOTAL_READS=$(grep -i "total.*reads" "$CHROMAP_LOG" | head -1 | grep -oE '[0-9,]+' | tr -d ',' || echo "N/A")
    echo "Total_input_reads\t$TOTAL_READS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    # Extract mapped reads
    MAPPED_READS=$(grep -i "mapped.*reads\|aligned.*reads" "$CHROMAP_LOG" | head -1 | grep -oE '[0-9,]+' | tr -d ',' || echo "N/A")
    echo "Mapped_reads\t$MAPPED_READS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    # Calculate mapping rate
    if [[ "$TOTAL_READS" != "N/A" && "$MAPPED_READS" != "N/A" ]]; then
        MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc)
        echo "Mapping_rate_percent\t$MAPPING_RATE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    else
        echo "Mapping_rate_percent\tN/A" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    fi
    
    # Extract duplicate rate if available
    DUPLICATE_RATE=$(grep -i "duplicate" "$CHROMAP_LOG" | grep -oE '[0-9.]+%' | head -1 || echo "N/A")
    echo "Duplicate_rate\t$DUPLICATE_RATE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt"
    
    echo "DEBUG: Chromap statistics extracted"
else
    echo "WARNING: Chromap log file not found: $CHROMAP_LOG"
fi

# =============================================================================
# FRAGMENT STATISTICS ANALYSIS
# PURPOSE: Analyze fragment distribution and barcode complexity metrics
# 
# KEY METRICS:
# - Total fragments: Overall fragment count from alignment
# - Unique barcodes: Number of distinct cell barcodes detected
# - Fragments per barcode: Distribution of fragments across barcodes
# 
# ANALYSIS STRATEGY:
# - Extract barcode information from fragment file (column 4)
# - Calculate barcode complexity and fragment distribution
# - Generate per-barcode statistics for cell quality assessment
# 
# QUALITY INDICATORS:
# - High unique barcode count suggests good cell capture
# - Balanced fragment distribution indicates uniform library preparation
# - These metrics help identify high-quality cells for downstream analysis
# 
# OUTPUT FILES:
# - fragment_summary.txt: Basic fragment and barcode counts
# - fragments_per_barcode_raw.txt: Per-barcode fragment counts for filtering
# =============================================================================
echo "DEBUG: Analyzing fragment statistics..."

# Basic fragment counts
TOTAL_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | wc -l)
UNIQUE_BARCODES=$(zcat "$FRAGMENTS_FILE" | cut -f4 | sort -u | wc -l)

echo "Total_fragments\t$TOTAL_FRAGMENTS" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt"
echo "Unique_barcodes\t$UNIQUE_BARCODES" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt"

# =============================================================================
# PER-BARCODE FRAGMENT DISTRIBUTION ANALYSIS
# PURPOSE: Calculate fragment distribution across individual cell barcodes
# 
# ANALYSIS WORKFLOW:
# 1. Extract barcodes from fragment file (column 4)
# 2. Count fragments per barcode using sort and uniq
# 3. Generate fragment count distribution for statistical analysis
# 4. Output raw fragment counts for downstream filtering
# 
# STATISTICAL APPLICATIONS:
# - Cell filtering based on fragment count thresholds
# - Quality assessment of single-cell capture efficiency
# - Identification of high-quality cells for peak calling
# - Distribution analysis for barcode complexity metrics
# =============================================================================
echo "DEBUG: Calculating per-barcode fragment statistics..."
zcat "$FRAGMENTS_FILE" | \
    cut -f4 | \
    sort | \
    uniq -c | \
    awk '{print $1}' | \
    sort -n > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt"

# Calculate summary statistics
FRAGMENT_STATS="$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_distribution.txt"
awk 'BEGIN{sum=0; count=0; min=""; max=""} 
     {
         sum+=$1; count++; 
         if(min=="" || $1<min) min=$1; 
         if(max=="" || $1>max) max=$1;
         values[count]=$1
     } 
     END{
         mean=sum/count
         asort(values)
         median=(count%2==1) ? values[int(count/2)+1] : (values[int(count/2)] + values[int(count/2)+1])/2
         print "Mean_fragments_per_barcode\t"mean
         print "Median_fragments_per_barcode\t"median
         print "Min_fragments_per_barcode\t"min
         print "Max_fragments_per_barcode\t"max
         print "Total_barcodes_with_fragments\t"count
     }' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" > "$FRAGMENT_STATS"

echo "DEBUG: Fragment distribution analysis completed"

# =============================================================================
# FRAGMENT LENGTH DISTRIBUTION ANALYSIS
# PURPOSE: Analyze fragment size distribution for chromatin accessibility assessment
# 
# FRAGMENT LENGTH SIGNIFICANCE:
# - Nucleosome-free regions: <147bp (accessible chromatin)
# - Mono-nucleosome: 147-294bp (single nucleosome)
# - Di-nucleosome: 294-441bp (two nucleosomes)
# - Larger fragments: Multi-nucleosome complexes
# 
# ANALYSIS WORKFLOW:
# 1. Calculate fragment lengths (end - start coordinates)
# 2. Sort lengths numerically for distribution analysis
# 3. Count frequency of each fragment length
# 4. Generate detailed length distribution table
# 
# QUALITY INDICATORS:
# - High proportion of nucleosome-free fragments indicates good accessibility
# - Clear nucleosome periodicity suggests proper ATAC-seq library preparation
# - Fragment length distribution informs peak calling strategy
# 
# OUTPUT FILES:
# - fragment_lengths_detailed.txt: Length frequency distribution table
# =============================================================================
echo "DEBUG: Analyzing fragment length distribution..."
zcat "$FRAGMENTS_FILE" | \
    awk '{print $3-$2}' | \
    sort -n | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_lengths_detailed.txt"

# Calculate fragment length statistics
FRAGMENT_LENGTH_STATS="$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_length_stats.txt"
zcat "$FRAGMENTS_FILE" | \
    awk '{len=$3-$2; sum+=len; sumsq+=len*len; lengths[NR]=len} 
         END{
             mean=sum/NR; 
             for(i=1;i<=NR;i++) variance+=(lengths[i]-mean)^2
             variance=variance/(NR-1)
             sd=sqrt(variance)
             asort(lengths)
             median=(NR%2==1) ? lengths[int(NR/2)+1] : (lengths[int(NR/2)] + lengths[int(NR/2)+1])/2
             q25=lengths[int(NR*0.25)]
             q75=lengths[int(NR*0.75)]
             
             print "Mean_fragment_length\t"mean
             print "Median_fragment_length\t"median
             print "SD_fragment_length\t"sd
             print "Q25_fragment_length\t"q25
             print "Q75_fragment_length\t"q75
             
             # Count fragments in key size ranges
             mono_nucleosome=0; di_nucleosome=0; long_fragments=0
             for(i=1;i<=NR;i++) {
                 if(lengths[i]<=150) mono_nucleosome++
                 else if(lengths[i]<=300) di_nucleosome++
                 else long_fragments++
             }
             
             print "Mono_nucleosome_fragments_<=150bp\t"mono_nucleosome
             print "Di_nucleosome_fragments_150-300bp\t"di_nucleosome
             print "Long_fragments_>300bp\t"long_fragments
             print "Mono_nucleosome_percentage\t"(mono_nucleosome*100/NR)
             print "Di_nucleosome_percentage\t"(di_nucleosome*100/NR)
         }' > "$FRAGMENT_LENGTH_STATS"

echo "DEBUG: Fragment length analysis completed"

# =============================================================================
# CHROMOSOME DISTRIBUTION ANALYSIS
# PURPOSE: Analyze fragment distribution across chromosomes and assess data quality
# 
# KEY METRICS:
# - Per-chromosome fragment counts and percentages
# - Mitochondrial DNA percentage (quality indicator)
# - Chromosome distribution balance assessment
# 
# ANALYSIS WORKFLOW:
# 1. Extract chromosome information from fragment file (column 1)
# 2. Count fragments per chromosome using sort and uniq
# 3. Calculate percentage distribution across chromosomes
# 4. Identify mitochondrial contamination levels
# 
# QUALITY INDICATORS:
# - Low mitochondrial percentage (<10%) indicates good nuclear enrichment
# - Balanced autosomal distribution suggests unbiased library preparation
# - High mitochondrial content may indicate cell stress or poor preparation
# 
# OUTPUT FILES:
# - chromosome_distribution.txt: Raw chromosome fragment counts
# - mito_stats.txt: Mitochondrial contamination metrics
# =============================================================================
echo "DEBUG: Analyzing chromosome distribution..."
zcat "$FRAGMENTS_FILE" | \
    cut -f1 | \
    sort | \
    uniq -c | \
    awk '{print $2"\t"$1}' | \
    sort -k1,1V > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_chromosome_distribution.txt"

# =============================================================================
# MITOCHONDRIAL CONTAMINATION ASSESSMENT
# PURPOSE: Calculate mitochondrial DNA percentage as quality control metric
# 
# MITOCHONDRIAL SIGNIFICANCE:
# - High mitochondrial percentage (>20%) indicates potential issues:
#   * Cell stress or damage during preparation
#   * Poor nuclear isolation
#   * Suboptimal ATAC-seq protocol execution
# - Low mitochondrial percentage (<10%) suggests good quality preparation
# 
# DETECTION STRATEGY:
# - Search for common mitochondrial chromosome names (chrM, MT)
# - Calculate percentage relative to total fragment count
# - Use flexible pattern matching for different genome annotations
# =============================================================================
MITO_FRAGMENTS=$(grep -E "^chrM|^MT" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_chromosome_distribution.txt" | awk '{sum+=$2} END{print sum+0}')
MITO_PERCENTAGE=$(echo "scale=2; $MITO_FRAGMENTS * 100 / $TOTAL_FRAGMENTS" | bc)

echo "Mitochondrial_fragments\t$MITO_FRAGMENTS" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt"
echo "Mitochondrial_percentage\t$MITO_PERCENTAGE" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt"

echo "DEBUG: Mitochondrial percentage: $MITO_PERCENTAGE%"

# =============================================================================
# READ QUALITY ASSESSMENT
# PURPOSE: Analyze read-level quality metrics and PCR duplication rates
# 
# KEY METRICS:
# - Total reads: Overall read count in BED format
# - Unique positions: Number of distinct genomic positions
# - PCR duplicate rate: Percentage of reads at identical positions
# 
# ANALYSIS STRATEGY:
# - Extract genomic coordinates (chr, start, end) from BED file
# - Count total reads and unique genomic positions
# - Calculate PCR duplication as reads sharing identical coordinates
# 
# QUALITY INDICATORS:
# - Low PCR duplicate rate (<30%) indicates good library complexity
# - High unique position count suggests diverse fragment coverage
# - These metrics complement fragment-based quality assessment
# 
# PCR DUPLICATE SIGNIFICANCE:
# - High duplication may indicate:
#   * Over-amplification during library preparation
#   * Low input cell numbers
#   * Poor library complexity
# - Affects downstream peak calling sensitivity and specificity
# =============================================================================
echo "DEBUG: Analyzing read quality from BED file..."
if [[ -f "$READS_FILE" ]]; then
    TOTAL_READS_BED=$(zcat "$READS_FILE" | wc -l)
    UNIQUE_READ_POSITIONS=$(zcat "$READS_FILE" | cut -f1-3 | sort -u | wc -l)
    
    echo "Total_reads_in_bed\t$TOTAL_READS_BED" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    echo "Unique_read_positions\t$UNIQUE_READ_POSITIONS" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    
    # Calculate PCR duplicate rate
    if [[ $TOTAL_READS_BED -gt 0 ]]; then
        DUPLICATE_RATE_CALC=$(echo "scale=4; (1 - $UNIQUE_READ_POSITIONS / $TOTAL_READS_BED) * 100" | bc)
        echo "Estimated_PCR_duplicate_rate_percent\t$DUPLICATE_RATE_CALC" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt"
    fi
fi

# =============================================================================
# BARCODE QUALITY ASSESSMENT
# PURPOSE: Evaluate cell barcode quality and classify cells by fragment content
# 
# QUALITY CATEGORIES:
# - High quality (≥1000 fragments): Likely high-quality cells
# - Medium quality (500-999 fragments): Moderate-quality cells
# - Low quality (<500 fragments): Low-quality cells or empty droplets
# 
# ASSESSMENT STRATEGY:
# - Use fragment count per barcode as primary quality metric
# - Classify barcodes into quality tiers for downstream filtering
# - Calculate distribution percentages for quality assessment
# 
# BIOLOGICAL INTERPRETATION:
# - High-quality barcodes represent viable cells with good accessibility
# - Medium-quality barcodes may represent stressed or partially lysed cells
# - Low-quality barcodes likely represent empty droplets or debris
# 
# FILTERING IMPLICATIONS:
# - High proportion of high-quality barcodes indicates successful experiment
# - Quality distribution informs cell filtering thresholds for peak calling
# - These metrics guide downstream analysis parameter selection
# =============================================================================
echo "DEBUG: Assessing barcode quality..."
# Barcode complexity
BARCODE_COMPLEXITY=$(echo "scale=4; $UNIQUE_BARCODES / $TOTAL_FRAGMENTS" | bc)
echo "Barcode_complexity\t$BARCODE_COMPLEXITY" > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"

# =============================================================================
# BARCODE COMPLEXITY CLASSIFICATION
# PURPOSE: Classify cell barcodes by fragment count into quality tiers
# 
# CLASSIFICATION THRESHOLDS:
# - High quality: ≥1000 fragments (robust cells for analysis)
# - Medium quality: 500-999 fragments (borderline cells)
# - Low quality: <500 fragments (likely empty droplets)
# 
# PERCENTAGE CALCULATIONS:
# - Calculate proportion of each quality tier
# - Assess overall experiment success rate
# - Inform cell filtering decisions for downstream analysis
# =============================================================================
BARCODES_GT_100=$(awk '$1>=100' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)
BARCODES_GT_500=$(awk '$1>=500' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)
BARCODES_GT_1000=$(awk '$1>=1000' "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragments_per_barcode_raw.txt" | wc -l)

echo "Barcodes_with_100plus_fragments\t$BARCODES_GT_100" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"
echo "Barcodes_with_500plus_fragments\t$BARCODES_GT_500" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"
echo "Barcodes_with_1000plus_fragments\t$BARCODES_GT_1000" >> "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt"

# =============================================================================
# COMPREHENSIVE QC SUMMARY GENERATION
# PURPOSE: Create a unified quality control report with all key metrics
# 
# SUMMARY COMPONENTS:
# - Alignment statistics: Mapping rates, duplicate rates
# - Fragment analysis: Counts, distributions, length characteristics
# - Chromosome distribution: Including mitochondrial contamination
# - Read quality: PCR duplication and position diversity
# - Barcode quality: Cell classification and quality distribution
# 
# REPORT STRUCTURE:
# - Organized sections for easy interpretation
# - Key metrics highlighted for quick assessment
# - Timestamp for tracking analysis runs
# 
# USAGE:
# - Primary QC document for experiment assessment
# - Input for downstream analysis parameter decisions
# - Reference for troubleshooting and optimization
# - Documentation for publication and sharing
# =============================================================================
echo "DEBUG: Creating comprehensive QC summary..."
cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_qc_summary.txt" << EOF
# Alignment Quality Control Summary for $SAMPLE
# Generated: $(date)
# Analysis performed before peak calling

=== ALIGNMENT STATISTICS ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" 2>/dev/null | grep -v "^#" || echo "Alignment_stats\tNot_available")

=== FRAGMENT SUMMARY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_fragment_summary.txt")

=== FRAGMENT DISTRIBUTION ===
$(cat "$FRAGMENT_STATS")

=== FRAGMENT LENGTH STATISTICS ===
$(cat "$FRAGMENT_LENGTH_STATS")

=== MITOCHONDRIAL CONTENT ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt")

=== READ QUALITY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_read_stats.txt" 2>/dev/null || echo "Read_stats\tNot_available")

=== BARCODE QUALITY ===
$(cat "$OUTPUT_DIR/alignment_qc/${SAMPLE}_barcode_quality.txt")
EOF

# =============================================================================
# R SCRIPT GENERATION FOR QC VISUALIZATION
# PURPOSE: Create comprehensive R script for generating quality control plots
# 
# PLOT TYPES:
# 1. Fragment length distribution: Shows nucleosome periodicity and accessibility
# 2. Fragments per barcode: Displays cell quality distribution
# 3. Chromosome distribution: Reveals mapping bias and mitochondrial content
# 
# VISUALIZATION FEATURES:
# - Reference lines for biological thresholds (nucleosome sizes, quality cutoffs)
# - Log-scale transformations for better data visualization
# - High-resolution output suitable for publication
# - Consistent styling and color schemes
# 
# SCRIPT CAPABILITIES:
# - Automatic data loading from QC analysis outputs
# - Error handling and directory creation
# - Flexible sample naming and output organization
# - Publication-ready plot generation
# 
# USAGE:
# - Run automatically if R is available
# - Can be executed manually for custom visualization
# - Generates PNG files in dedicated plots directory
# =============================================================================
echo "DEBUG: Creating R plotting script..."
cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_plot_alignment_qc.R" << 'EOF'
#!/usr/bin/env Rscript
# =============================================================================
# QC VISUALIZATION SCRIPT
# PURPOSE: Generate comprehensive quality control plots for ATAC-seq alignment
# 
# REQUIRED LIBRARIES:
# - ggplot2: Advanced plotting and visualization
# - dplyr: Data manipulation and filtering
# - gridExtra: Multi-panel plot arrangement
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("Please provide sample name as argument")
}

sample_name <- args[1]

suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
})

setwd("alignment_qc")

# =============================================================================
# FRAGMENT LENGTH DISTRIBUTION PLOT
# PURPOSE: Visualize fragment size distribution and nucleosome periodicity
# 
# FEATURES:
# - Reference lines at nucleosome boundaries (147bp, 294bp)
# - Focus on 0-1000bp range for clarity
# - Annotations for biological interpretation
# =============================================================================
if(file.exists(paste0(sample_name, "_fragment_lengths_detailed.txt"))) {
    frag_lengths <- read.table(paste0(sample_name, "_fragment_lengths_detailed.txt"), 
                              col.names = c("Length", "Count"))
    
    p1 <- ggplot(frag_lengths %>% filter(Length <= 1000), aes(x = Length, y = Count)) +
        geom_line(color = "blue", alpha = 0.8) +
        geom_vline(xintercept = 147, color = "red", linetype = "dashed", alpha = 0.7) +
        geom_vline(xintercept = 294, color = "orange", linetype = "dashed", alpha = 0.7) +
        annotate("text", x = 147, y = max(frag_lengths$Count[frag_lengths$Length <= 1000]) * 0.8, 
                 label = "Mono-nucleosome", angle = 90, vjust = -0.5, size = 3) +
        annotate("text", x = 294, y = max(frag_lengths$Count[frag_lengths$Length <= 1000]) * 0.8, 
                 label = "Di-nucleosome", angle = 90, vjust = -0.5, size = 3) +
        labs(title = paste("Fragment Length Distribution -", sample_name),
             x = "Fragment Length (bp)", y = "Count") +
        theme_minimal()
    
    ggsave(paste0(sample_name, "_fragment_length_distribution.png"), p1, width = 10, height = 6, dpi = 300)
}

# =============================================================================
# FRAGMENTS PER BARCODE DISTRIBUTION PLOT
# PURPOSE: Visualize cell quality distribution and filtering thresholds
# 
# FEATURES:
# - Log10 scale for better visualization of wide range
# - Reference lines at quality thresholds (100, 500, 1000 fragments)
# - Histogram to show barcode count distribution
# =============================================================================
if(file.exists(paste0(sample_name, "_fragments_per_barcode_raw.txt"))) {
    frags_per_bc <- read.table(paste0(sample_name, "_fragments_per_barcode_raw.txt"), 
                              col.names = "Fragment_count")
    
    p2 <- ggplot(frags_per_bc, aes(x = Fragment_count)) +
        geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue") +
        scale_x_log10() +
        geom_vline(xintercept = 100, color = "red", linetype = "dashed") +
        geom_vline(xintercept = 500, color = "orange", linetype = "dashed") +
        geom_vline(xintercept = 1000, color = "green", linetype = "dashed") +
        labs(title = paste("Fragments per Barcode Distribution -", sample_name),
             x = "Fragments per Barcode (log10)", y = "Number of Barcodes") +
        theme_minimal()
    
    ggsave(paste0(sample_name, "_fragments_per_barcode.png"), p2, width = 10, height = 6, dpi = 300)
}

# =============================================================================
# CHROMOSOME DISTRIBUTION PLOT
# PURPOSE: Visualize fragment distribution across chromosomes
# 
# FEATURES:
# - Focus on main chromosomes for clarity
# - Ordered chromosome display for easy interpretation
# - Bar chart for clear quantitative comparison
# =============================================================================
if(file.exists(paste0(sample_name, "_chromosome_distribution.txt"))) {
    chrom_dist <- read.table(paste0(sample_name, "_chromosome_distribution.txt"), 
                            col.names = c("Chromosome", "Count"))
    
    # Filter for main chromosomes and sort
    main_chroms <- chrom_dist %>% 
        filter(grepl("^chr[0-9XY]+$", Chromosome)) %>%
        mutate(Chromosome = factor(Chromosome, 
                                  levels = paste0("chr", c(1:19, "X", "Y"))))
    
    p3 <- ggplot(main_chroms, aes(x = Chromosome, y = Count)) +
        geom_bar(stat = "identity", alpha = 0.7, fill = "lightcoral") +
        labs(title = paste("Fragment Distribution by Chromosome -", sample_name),
             x = "Chromosome", y = "Fragment Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(sample_name, "_chromosome_distribution.png"), p3, width = 12, height = 6, dpi = 300)
}

cat("Alignment QC plots generated successfully!\n")
EOF

chmod +x "$OUTPUT_DIR/alignment_qc/${SAMPLE}_plot_alignment_qc.R"

# =============================================================================
# R SCRIPT EXECUTION FOR PLOT GENERATION
# PURPOSE: Automatically generate QC plots if R environment is available
# 
# EXECUTION STRATEGY:
# - Check for R availability in system PATH
# - Execute R script with proper arguments and logging
# - Provide fallback instructions if R is unavailable
# 
# ERROR HANDLING:
# - Capture R script output and errors in log file
# - Continue pipeline execution even if plotting fails
# - Provide clear instructions for manual plot generation
# 
# OUTPUT:
# - High-resolution PNG plots in dedicated directory
# - Execution log for troubleshooting
# - Clear status messages for user feedback
# 
# BENEFITS:
# - Immediate visual QC assessment
# - Publication-ready plot generation
# - Automated workflow integration
# =============================================================================
if command -v Rscript &> /dev/null; then
    echo "DEBUG: Generating alignment QC plots..."
    cd "$OUTPUT_DIR"
    Rscript "alignment_qc/${SAMPLE}_plot_alignment_qc.R" "$SAMPLE" || echo "WARNING: R plotting failed"
    cd - > /dev/null
else
    echo "WARNING: Rscript not available. Plots not generated."
fi

# =============================================================================
# QUALITY ASSESSMENT AND RECOMMENDATIONS GENERATION
# PURPOSE: Provide automated quality evaluation and actionable recommendations
# 
# ASSESSMENT CRITERIA:
# - Mapping rate: Alignment efficiency (>80% excellent, >70% good, >60% moderate)
# - Mono-nucleosome %: Chromatin accessibility (>40% excellent, >30% good, >20% moderate)
# - Mitochondrial %: Contamination level (<10% excellent, <20% good, <30% moderate)
# - High-quality barcodes: Cell recovery (>20% excellent, >15% good, >10% moderate)
# 
# RECOMMENDATION SYSTEM:
# - Automated quality scoring with visual indicators
# - Specific thresholds based on ATAC-seq best practices
# - Actionable next steps for different quality levels
# - Integration guidance for downstream analysis
# 
# OUTPUT FEATURES:
# - Clear quality indicators (✓ excellent/good, ⚠ moderate, ✗ poor)
# - Specific recommendations for improvement
# - Next step guidance for pipeline continuation
# - Troubleshooting suggestions for poor quality samples
# =============================================================================
echo "DEBUG: Generating quality assessment recommendations..."

# Read key metrics for recommendations
MAPPING_RATE_NUM=$(grep "Mapping_rate_percent" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_stats.txt" 2>/dev/null | cut -f2 || echo "0")
MONO_PERCENTAGE=$(grep "Mono_nucleosome_percentage" "$FRAGMENT_LENGTH_STATS" | cut -f2)
MITO_PERC=$(grep "Mitochondrial_percentage" "$OUTPUT_DIR/alignment_qc/${SAMPLE}_mito_stats.txt" | cut -f2)
HIGH_QUALITY_BARCODES=$BARCODES_GT_500

cat > "$OUTPUT_DIR/alignment_qc/${SAMPLE}_quality_recommendations.txt" << EOF
# Quality Assessment and Recommendations for $SAMPLE
# Generated: $(date)

=== OVERALL ASSESSMENT ===
Sample: $SAMPLE
Analysis Date: $(date)

=== KEY METRICS ===
Mapping Rate: ${MAPPING_RATE_NUM}%
Mono-nucleosome Percentage: ${MONO_PERCENTAGE}%
Mitochondrial Percentage: ${MITO_PERC}%
High-quality Barcodes (>500 frags): $HIGH_QUALITY_BARCODES

=== QUALITY RECOMMENDATIONS ===

MAPPING RATE:
$(if (( $(echo "$MAPPING_RATE_NUM > 60" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: Mapping rate > 60% indicates good alignment quality"
elif (( $(echo "$MAPPING_RATE_NUM > 40" | bc -l 2>/dev/null || echo 0) )); then
    echo "⚠️  FAIR: Mapping rate 40-60% is acceptable but could be improved"
else
    echo "❌ POOR: Mapping rate < 40% suggests alignment issues"
fi)

FRAGMENT LENGTH DISTRIBUTION:
$(if (( $(echo "$MONO_PERCENTAGE > 40" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: High mono-nucleosome percentage indicates good nucleosome structure"
else
    echo "⚠️  CHECK: Low mono-nucleosome percentage - verify fragment size selection"
fi)

MITOCHONDRIAL CONTENT:
$(if (( $(echo "$MITO_PERC < 20" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ GOOD: Low mitochondrial percentage indicates good nuclear enrichment"
elif (( $(echo "$MITO_PERC < 30" | bc -l 2>/dev/null || echo 0) )); then
    echo "⚠️  FAIR: Moderate mitochondrial content - still usable"
else
    echo "❌ POOR: High mitochondrial content suggests poor nuclear enrichment"
fi)

CELL QUALITY:
$(if [[ $HIGH_QUALITY_BARCODES -gt 1000 ]]; then
    echo "✅ GOOD: >1000 high-quality cells detected"
elif [[ $HIGH_QUALITY_BARCODES -gt 500 ]]; then
    echo "⚠️  FAIR: 500-1000 high-quality cells detected"
else
    echo "❌ POOR: <500 high-quality cells detected"
fi)

=== RECOMMENDATION FOR PEAK CALLING ===
$(
overall_score=0
if (( $(echo "$MAPPING_RATE_NUM > 50" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if (( $(echo "$MONO_PERCENTAGE > 35" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if (( $(echo "$MITO_PERC < 25" | bc -l 2>/dev/null || echo 0) )); then ((overall_score++)); fi
if [[ $HIGH_QUALITY_BARCODES -gt 500 ]]; then ((overall_score++)); fi

if [[ $overall_score -ge 3 ]]; then
    echo "✅ PROCEED: Data quality is sufficient for peak calling"
    echo "   Recommended to continue with step 5 (05_call_peaks.sh)"
elif [[ $overall_score -ge 2 ]]; then
    echo "⚠️  PROCEED WITH CAUTION: Data quality is marginal"
    echo "   Consider reviewing alignment parameters or filtering thresholds"
    echo "   You may still proceed with peak calling but monitor results carefully"
else
    echo "❌ REVIEW REQUIRED: Data quality issues detected"
    echo "   Recommend reviewing alignment parameters before peak calling"
    echo "   Consider re-processing with adjusted settings"
fi
)

=== FILES GENERATED ===
- Comprehensive QC summary: ${SAMPLE}_alignment_qc_summary.txt
- Fragment statistics: ${SAMPLE}_fragment_distribution.txt
- Length statistics: ${SAMPLE}_fragment_length_stats.txt
- Chromosome distribution: ${SAMPLE}_chromosome_distribution.txt
- Barcode quality: ${SAMPLE}_barcode_quality.txt
- QC plots: ${SAMPLE}_*.png (if R available)
- This recommendation file: ${SAMPLE}_quality_recommendations.txt
EOF

# =============================================================================
# PIPELINE COMPLETION AND SUMMARY
# PURPOSE: Provide comprehensive completion summary and next step guidance
# 
# COMPLETION SUMMARY:
# - Display key quality metrics for quick assessment
# - List all generated output files and their purposes
# - Provide clear next step instructions for pipeline continuation
# 
# KEY METRICS DISPLAY:
# - Total fragments: Overall data volume
# - Unique barcodes: Cell recovery estimate
# - Mapping rate: Alignment efficiency
# - Mitochondrial content: Contamination level
# - High-quality cells: Usable cell count
# 
# OUTPUT FILES SUMMARY:
# - QC summary: Comprehensive metrics overview
# - Quality recommendations: Automated assessment and guidance
# - Fragment analysis: Detailed distributions and statistics
# - Visualization: R scripts and plots for visual assessment
# 
# NEXT STEPS GUIDANCE:
# - Quality review workflow
# - Decision points for pipeline continuation
# - File locations for easy access
# =============================================================================
echo "Quality Assessment Summary:"
echo "=========================="
echo "Sample: $SAMPLE"
echo "Total Fragments: $(printf "%'d" $TOTAL_FRAGMENTS)"
echo "Unique Barcodes: $(printf "%'d" $UNIQUE_BARCODES)"
echo "Mapping Rate: ${MAPPING_RATE_NUM}%"
echo "Mitochondrial Content: ${MITO_PERC}%"
echo "High-Quality Cells (>500 frags): $(printf "%'d" $HIGH_QUALITY_BARCODES)"
echo ""
echo "📊 Detailed QC report: $OUTPUT_DIR/alignment_qc/${SAMPLE}_alignment_qc_summary.txt"
echo "📋 Quality recommendations: $OUTPUT_DIR/alignment_qc/${SAMPLE}_quality_recommendations.txt"

echo "Output files created in: $OUTPUT_DIR/alignment_qc/"
echo "  - Alignment QC summary: ${SAMPLE}_alignment_qc_summary.txt"
echo "  - Quality recommendations: ${SAMPLE}_quality_recommendations.txt"
echo "  - Fragment statistics: ${SAMPLE}_fragment_distribution.txt"
echo "  - Fragment length stats: ${SAMPLE}_fragment_length_stats.txt"
echo "  - Chromosome distribution: ${SAMPLE}_chromosome_distribution.txt"
echo "  - R plotting script: ${SAMPLE}_plot_alignment_qc.R"

echo "========================================="
echo "Step 4b complete for $SAMPLE"
echo "Alignment quality control analysis finished"
echo "End time: $(date)"
echo "========================================="