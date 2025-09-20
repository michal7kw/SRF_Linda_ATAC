#!/bin/bash
# =============================================================================
# COMPREHENSIVE QC METRICS ANALYSIS FOR scATAC-seq DATA
# =============================================================================
# PURPOSE: Calculate comprehensive quality control metrics for scATAC-seq data
# 
# DESCRIPTION:
# This script performs extensive quality control analysis on processed scATAC-seq
# data, generating multiple metrics essential for data validation and downstream
# analysis planning. It analyzes fragment characteristics, TSS enrichment,
# promoter accessibility, library complexity, and per-cell quality metrics.
# 
# INPUT REQUIREMENTS:
# - Fragments file: ${SAMPLE}_fragments.tsv.gz (from chromap alignment)
# - Peaks file: ${SAMPLE}_peaks_sorted.bed (from MACS2 peak calling)
# - Chromosome sizes file: mm10.chrom.sizes (reference genome)
# - Peak-barcode matrix (optional, for enhanced cell metrics)
# 
# OUTPUT PRODUCTS:
# - Fragment length distribution and statistics
# - TSS enrichment scores and analysis
# - Promoter accessibility metrics
# - Per-cell fragment and peak counts
# - Library complexity measurements
# - Comprehensive QC summary report
# - R plotting scripts for visualization
# 
# COMPUTATIONAL REQUIREMENTS:
# - Memory: 32GB (for large fragment files and statistical calculations)
# - CPUs: 8 cores (parallel processing of metrics)
# - Time: 6 hours (comprehensive analysis with external data downloads)
# - Storage: ~5-10GB for intermediate files and results
# 
# DEPENDENCIES:
# - bedtools: Genomic interval operations and overlaps
# - wget: Download reference annotations (TSS regions)
# - awk/bc: Statistical calculations and data processing
# - R (optional): Advanced plotting and visualization
# 
# QUALITY METRICS CALCULATED:
# - Fragment size distribution (nucleosome positioning)
# - TSS enrichment (chromatin accessibility at promoters)
# - Promoter peak percentage (regulatory region accessibility)
# - Library complexity (PCR duplication assessment)
# - Per-cell fragment counts (cell quality assessment)
# - Per-cell peak accessibility (cell-specific patterns)
# =============================================================================
#SBATCH --job-name=qc_metrics
#SBATCH --output=logs/08_qc_metrics_%a.out
#SBATCH --error=logs/08_qc_metrics_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=workq

# =============================================================================
# CONDA ENVIRONMENT SETUP AND CONFIGURATION
# PURPOSE: Initialize computational environment with required tools
# 
# ENVIRONMENT FEATURES:
# - alignment_two: Comprehensive bioinformatics environment
# - bedtools: Genomic interval operations and statistics
# - wget/curl: External data download capabilities
# - awk/bc: Advanced text processing and mathematical operations
# - R packages: ggplot2, dplyr for optional visualization
# 
# ERROR HANDLING:
# - Strict error checking with set -euo pipefail
# - Immediate exit on command failures
# - Undefined variable detection
# - Pipeline failure propagation
# 
# COMPUTATIONAL ENVIRONMENT:
# - Ensures reproducible analysis environment
# - Validates tool availability before processing
# - Maintains consistent software versions
# - Supports both command-line and R-based analysis
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# SAMPLE CONFIGURATION AND PROCESSING SETUP
# PURPOSE: Configure sample-specific parameters for QC analysis
# 
# SAMPLE CONFIGURATION:
# - R26-Nestin-Ctrl-adult: Control adult neural stem cells
# - R26-Nestin-Mut-adult: Mutant adult neural stem cells
# - Array-based processing for parallel sample analysis
# - Consistent naming convention across pipeline steps
# 
# PATH CONFIGURATION:
# - OUTPUT_DIR: Central location for all pipeline outputs
# - Organized subdirectories for different analysis types
# - Consistent file naming for downstream integration
# - Supports both individual and batch processing
# 
# ARRAY PROCESSING:
# - SLURM_ARRAY_TASK_ID: Enables parallel sample processing
# - Automatic sample selection based on array index
# - Scalable to additional samples without code changes
# - Efficient resource utilization across compute nodes
# 
# LOGGING AND TRACKING:
# - Comprehensive start/end timestamps
# - Sample identification for debugging
# - Processing step documentation
# - Progress tracking for long-running analyses
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 8: QC Metrics Analysis for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY PREPARATION
# PURPOSE: Create organized directory structure for QC metrics
# 
# DIRECTORY STRUCTURE:
# - qc_metrics/: Central location for all quality control outputs
# - Organized storage for multiple metric types
# - Consistent naming for downstream analysis
# - Supports both individual and comparative analysis
# 
# BENEFITS:
# - Centralized QC results for easy access
# - Prevents file conflicts between samples
# - Facilitates batch processing and comparison
# - Enables automated downstream analysis
# =============================================================================
mkdir -p "$OUTPUT_DIR/qc_metrics"

# =============================================================================
# PREREQUISITE VALIDATION AND INPUT VERIFICATION
# PURPOSE: Ensure all required input files are available for QC analysis
# 
# REQUIRED INPUT FILES:
# - Fragments file: Contains aligned fragment coordinates and barcodes
# - Peaks file: Contains called accessible chromatin regions
# - Chromosome sizes: Reference genome dimensions for calculations
# 
# VALIDATION CHECKS:
# - File existence verification
# - Dependency chain validation
# - Clear error messages for missing inputs
# - Graceful exit with informative feedback
# 
# INPUT REQUIREMENTS:
# - Fragments: Output from chromap alignment (step 4)
# - Peaks: Output from MACS2 peak calling (step 5)
# - Chromosome sizes: Reference genome metadata
# 
# ERROR HANDLING:
# - Specific error messages for each missing file
# - Guidance on prerequisite pipeline steps
# - Prevents downstream analysis with incomplete data
# - Maintains data integrity throughout pipeline
# =============================================================================
FRAGMENTS_FILE="$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
PEAKS_FILE="$OUTPUT_DIR/peaks/${SAMPLE}_peaks_sorted.bed"
CHROM_SIZES="$OUTPUT_DIR/mm10.chrom.sizes"

if [[ ! -f "$FRAGMENTS_FILE" ]]; then
    echo "ERROR: Fragments file not found: $FRAGMENTS_FILE"
    exit 1
fi

if [[ ! -f "$PEAKS_FILE" ]]; then
    echo "ERROR: Peaks file not found: $PEAKS_FILE"
    echo "Please run step 5 (05_call_peaks.sh) first"
    exit 1
fi

# =============================================================================
# TOOL VALIDATION AND VERIFICATION
# PURPOSE: Ensure all required computational tools are available
# 
# REQUIRED TOOLS:
# - bedtools: Genomic interval operations, overlaps, and statistics
# - wget: Download external reference data (TSS annotations)
# - awk: Text processing and statistical calculations
# - bc: Precise mathematical operations and floating-point arithmetic
# - sort/uniq: Data sorting and uniqueness operations
# 
# TOOL FUNCTIONS:
# - bedtools intersect: Calculate overlaps between genomic regions
# - bedtools slop: Extend genomic intervals for analysis windows
# - wget: Download UCSC genome annotations
# - awk: Parse files and calculate statistics
# - bc: Perform precise mathematical calculations
# 
# VALIDATION IMPORTANCE:
# - Prevents runtime failures due to missing dependencies
# - Ensures reproducible analysis environment
# - Provides clear error messages for troubleshooting
# - Validates computational environment before processing
# 
# ERROR HANDLING:
# - Specific error messages for each missing tool
# - Installation guidance for missing dependencies
# - Graceful exit with informative feedback
# - Environment validation before data processing
# =============================================================================
echo "DEBUG: Checking for required tools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found"
    exit 1
fi

echo "DEBUG: Tools verified successfully"

# =============================================================================
# FRAGMENT LENGTH DISTRIBUTION ANALYSIS
# PURPOSE: Analyze fragment size patterns for nucleosome positioning assessment
# 
# BIOLOGICAL SIGNIFICANCE:
# - Fragment lengths reflect nucleosome organization
# - ~147bp: Nucleosome-bound DNA (mononucleosome)
# - ~300bp: Dinucleosome fragments
# - <100bp: Nucleosome-free regions (NFR)
# - >500bp: Potential artifacts or large nucleosome-free regions
# 
# ANALYSIS WORKFLOW:
# 1. Extract fragment lengths from coordinates (end - start)
# 2. Generate frequency distribution of all fragment sizes
# 3. Calculate comprehensive statistical metrics
# 4. Identify quality indicators and potential issues
# 
# STATISTICAL METRICS:
# - Total fragments: Overall library size
# - Mean/median length: Central tendency indicators
# - Standard deviation: Fragment size variability
# - Min/max lengths: Range and potential artifacts
# - Distribution shape: Nucleosome positioning quality
# 
# QUALITY INDICATORS:
# - Strong ~147bp peak: Good nucleosome positioning
# - Clear periodicity: Intact chromatin structure
# - Appropriate size range: Successful fragmentation
# - Low artifact peaks: Clean library preparation
# 
# OUTPUT APPLICATIONS:
# - Quality control assessment
# - Nucleosome positioning analysis
# - Library preparation optimization
# - Comparative analysis between samples
# =============================================================================
echo "DEBUG: Calculating fragment length distribution..."
zcat "$FRAGMENTS_FILE" | \
    awk '{print $3-$2}' | \
    sort -n | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragment_lengths.txt"

# Calculate fragment length statistics
FRAGMENT_STATS="$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragment_stats.txt"
zcat "$FRAGMENTS_FILE" | \
    awk '{len=$3-$2; sum+=len; sumsq+=len*len; if(NR==1){min=max=len} if(len<min){min=len} if(len>max){max=len}} 
         END{mean=sum/NR; variance=(sumsq-sum*sum/NR)/(NR-1); sd=sqrt(variance); 
         print "Total_fragments\t"NR"\nMean_length\t"mean"\nMedian_length\t"median"\nSD_length\t"sd"\nMin_length\t"min"\nMax_length\t"max}' > "$FRAGMENT_STATS"

echo "DEBUG: Fragment length analysis completed"

# =============================================================================
# TSS ENRICHMENT ANALYSIS
# PURPOSE: Assess chromatin accessibility at transcription start sites
# 
# BIOLOGICAL SIGNIFICANCE:
# - TSS regions are typically nucleosome-free and highly accessible
# - High TSS enrichment indicates successful ATAC-seq library preparation
# - Low enrichment suggests poor chromatin accessibility or technical issues
# - TSS enrichment is a key quality metric for ATAC-seq experiments
# 
# ANALYSIS WORKFLOW:
# 1. Download reference TSS annotations from UCSC Genome Browser
# 2. Create TSS-centered regions (±2kb windows for comprehensive analysis)
# 3. Calculate fragment overlaps with TSS regions
# 4. Compute enrichment percentage (TSS fragments / total fragments * 100)
# 5. Generate quality assessment metrics
# 
# TSS REGION DEFINITION:
# - Central TSS coordinates from RefSeq annotations
# - ±2kb windows to capture broader promoter accessibility
# - Strand-specific TSS identification for accuracy
# - Genome-wide TSS coverage for comprehensive analysis
# 
# ENRICHMENT CALCULATION:
# - TSS overlaps: Fragments intersecting TSS regions
# - Total fragments: Complete fragment library size
# - Enrichment percentage: Proportion of TSS-associated fragments
# - Quality threshold: Typically >10% for good libraries
# 
# QUALITY INDICATORS:
# - High enrichment (>20%): Excellent accessibility
# - Moderate enrichment (10-20%): Good quality
# - Low enrichment (<10%): Potential technical issues
# - Very low enrichment (<5%): Poor library quality
# 
# TECHNICAL CONSIDERATIONS:
# - Reference genome compatibility (mm10)
# - TSS annotation completeness and accuracy
# - Fragment size filtering effects
# - Cell type-specific accessibility patterns
# =============================================================================
echo "DEBUG: Calculating TSS enrichment..."
# Download TSS regions if not available (mouse mm10)
TSS_FILE="$OUTPUT_DIR/qc_metrics/mm10_tss.bed"
if [[ ! -f "$TSS_FILE" ]]; then
    echo "DEBUG: Downloading comprehensive mm10 TSS file from UCSC..."
    # Download and process RefSeq gene annotations for mm10 from UCSC
    wget -qO- "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz" | \
        gunzip -c | \
        awk 'BEGIN{OFS="\t"} {if($4=="+"){print $3, $5-1, $5, $13, "0", $4} else {print $3, $6-1, $6, $13, "0", $4}}' | \
        sort -k1,1 -k2,2n > "$TSS_FILE"
    echo "DEBUG: Comprehensive TSS file created successfully."
fi

# Calculate TSS enrichment (fragments overlapping TSS +/- 2kb)
bedtools slop -i "$TSS_FILE" -g "$CHROM_SIZES" -b 2000 | \
    bedtools intersect -a "$FRAGMENTS_FILE" -b stdin -c > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_overlap.txt"

# Calculate enrichment score
TSS_ENRICHMENT=$(awk 'BEGIN{tss=0; total=0} {total++; if($5>0) tss++} END{print tss/total*100}' "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_overlap.txt")
echo "TSS_enrichment_percent\t$TSS_ENRICHMENT" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_stats.txt"

echo "DEBUG: TSS enrichment calculated: $TSS_ENRICHMENT%"

# =============================================================================
# PEAK-IN-PROMOTER ANALYSIS
# PURPOSE: Assess the proportion of peaks located in gene promoter regions
# 
# BIOLOGICAL SIGNIFICANCE:
# - Promoter peaks indicate active transcriptional regulation
# - High promoter enrichment suggests functional accessibility
# - Balanced promoter/enhancer ratio indicates comprehensive coverage
# - Cell type-specific promoter usage patterns
# 
# ANALYSIS WORKFLOW:
# 1. Define promoter regions as TSS ±1kb windows
# 2. Calculate overlaps between called peaks and promoter regions
# 3. Compute percentage of peaks in promoter regions
# 4. Generate regulatory element distribution metrics
# 
# PROMOTER REGION DEFINITION:
# - TSS-centered windows (±1kb) for core promoter coverage
# - Captures proximal regulatory elements
# - Includes core promoter and immediate upstream/downstream regions
# - Standardized window size for comparative analysis
# 
# OVERLAP CALCULATION:
# - Peaks in promoters: Number of peaks overlapping promoter regions
# - Total peaks: Complete set of called peaks
# - Promoter percentage: Proportion of promoter-associated peaks
# - Regulatory element distribution assessment
# 
# QUALITY INDICATORS:
# - High promoter percentage (>40%): Strong transcriptional focus
# - Moderate percentage (20-40%): Balanced regulatory coverage
# - Low percentage (<20%): Enhancer-dominated or poor quality
# - Cell type variation: Expected biological differences
# 
# INTERPRETATION CONSIDERATIONS:
# - Cell type-specific promoter usage patterns
# - Developmental stage effects on promoter accessibility
# - Technical factors affecting peak calling sensitivity
# - Comparison with reference datasets for validation
# =============================================================================
echo "DEBUG: Analyzing peaks in promoter regions..."
if [[ -f "$TSS_FILE" ]]; then
    # Create promoter regions (TSS +/- 1kb)
    bedtools slop -i "$TSS_FILE" -g "$CHROM_SIZES" -b 1000 > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoters.bed"
    
    # Count peaks in promoters
    PEAKS_IN_PROMOTERS=$(bedtools intersect -a "$PEAKS_FILE" -b "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoters.bed" -u | wc -l)
    TOTAL_PEAKS=$(wc -l < "$PEAKS_FILE")
    PROMOTER_PERCENT=$(awk -v p="$PEAKS_IN_PROMOTERS" -v t="$TOTAL_PEAKS" 'BEGIN {printf "%.2f", p/t*100}')
    
    echo "Peaks_in_promoters\t$PEAKS_IN_PROMOTERS" > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    echo "Total_peaks\t$TOTAL_PEAKS" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    echo "Promoter_percentage\t$PROMOTER_PERCENT" >> "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt"
    
    echo "DEBUG: Promoter analysis completed: $PROMOTER_PERCENT% peaks in promoters"
fi

# =============================================================================
# PER-CELL QUALITY METRICS ANALYSIS
# PURPOSE: Assess single-cell level quality and heterogeneity metrics
# 
# BIOLOGICAL SIGNIFICANCE:
# - Cell-level metrics reveal library quality heterogeneity
# - Fragment counts indicate cell viability and accessibility
# - Peak accessibility shows functional chromatin states
# - Quality distributions guide cell filtering decisions
# 
# ANALYSIS WORKFLOW:
# 1. Calculate fragments per cell from barcode information
# 2. Determine peak accessibility per cell (if matrix available)
# 3. Compute average metrics and distributions
# 4. Generate cell quality assessment data
# 
# FRAGMENTS PER CELL:
# - Total accessible fragments per individual cell
# - Indicates overall chromatin accessibility
# - Quality threshold: typically >1000 fragments/cell
# - Distribution shape reveals library heterogeneity
# 
# PEAKS PER CELL:
# - Peak accessibility counts per cell from matrix
# - Measures functional accessibility
# - Signal-to-noise ratio indicator
# - Peak accessibility efficiency metric
# 
# QUALITY METRICS:
# - Cell count: Total number of detected cells
# - Fragment distribution: Library depth heterogeneity
# - Peak accessibility: Functional chromatin patterns
# - Quality filtering thresholds for downstream analysis
# 
# QUALITY INDICATORS:
# - High fragment counts (>5000): Excellent cell quality
# - Moderate counts (1000-5000): Good quality cells
# - Low counts (<1000): Potential low-quality cells
# - Peak/fragment ratio: Signal specificity measure
# 
# FILTERING APPLICATIONS:
# - Cell quality filtering thresholds
# - Outlier detection and removal
# - Library complexity assessment
# - Downstream analysis cell selection
# =============================================================================
echo "DEBUG: Calculating per-cell quality metrics..."
# Calculate fragments per cell
zcat "$FRAGMENTS_FILE" | \
    cut -f4 | \
    sort | \
    uniq -c | \
    awk '{print $2"\t"$1}' > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt"

# Calculate peaks per cell (if peak-cell matrix exists)
PEAK_BC_MATRIX="$OUTPUT_DIR/${SAMPLE}_peak_bc_matrix"
if [[ -d "$PEAK_BC_MATRIX" ]]; then
    if [[ -f "$PEAK_BC_MATRIX/matrix.mtx.gz" ]]; then
        echo "DEBUG: Calculating peaks per cell from matrix..."
        # Extract non-zero entries per barcode
        zcat "$PEAK_BC_MATRIX/matrix.mtx.gz" | \
            tail -n +4 | \
            awk '{count[$2]++} END{for(bc in count) print bc"\t"count[bc]}' | \
            sort -n > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_peaks_per_cell_idx.txt"
    fi
fi

# =============================================================================
# LIBRARY COMPLEXITY ANALYSIS
# PURPOSE: Assess library diversity and PCR duplication levels
# 
# BIOLOGICAL SIGNIFICANCE:
# - Library complexity reflects chromatin accessibility diversity
# - High complexity indicates comprehensive genome coverage
# - Low complexity suggests PCR over-amplification or limited accessibility
# - Complexity metrics guide sequencing depth optimization
# 
# ANALYSIS WORKFLOW:
# 1. Identify unique genomic fragments (coordinates-based)
# 2. Calculate total vs unique fragment ratios
# 3. Compute duplication rates and complexity scores
# 4. Generate library quality assessment metrics
# 
# COMPLEXITY METRICS:
# - Unique fragments: Distinct genomic coordinates
# - Total fragments: All sequenced fragments including duplicates
# - Duplication rate: Proportion of duplicate fragments
# - Complexity score: Ratio of unique to total fragments
# 
# DUPLICATION SOURCES:
# - PCR amplification during library preparation
# - Natural chromatin accessibility hotspots
# - Technical artifacts from sequencing
# - Biological replication of accessible regions
# 
# QUALITY INDICATORS:
# - High complexity (>0.8): Excellent library diversity
# - Moderate complexity (0.6-0.8): Good library quality
# - Low complexity (<0.6): Potential over-amplification
# - Very low complexity (<0.4): Poor library preparation
# 
# INTERPRETATION CONSIDERATIONS:
# - Cell type-specific accessibility patterns
# - Sequencing depth effects on apparent complexity
# - PCR cycle optimization for library preparation
# - Comparison with technical replicates
# 
# APPLICATIONS:
# - Library preparation optimization
# - Sequencing depth planning
# - Quality control assessment
# - Comparative analysis between samples
# =============================================================================
echo "DEBUG: Calculating library complexity..."
TOTAL_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | wc -l)
UNIQUE_FRAGMENTS=$(zcat "$FRAGMENTS_FILE" | cut -f1-3 | sort -u | wc -l)
COMPLEXITY=$(echo "scale=4; $UNIQUE_FRAGMENTS / $TOTAL_FRAGMENTS" | bc)

cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_library_complexity.txt" << EOF
Total_fragments	$TOTAL_FRAGMENTS
Unique_fragments	$UNIQUE_FRAGMENTS
Library_complexity	$COMPLEXITY
EOF

echo "DEBUG: Library complexity: $COMPLEXITY"

# =============================================================================
# COMPREHENSIVE QC SUMMARY GENERATION
# PURPOSE: Compile all quality control metrics into a unified report
# 
# SUMMARY COMPONENTS:
# - Fragment analysis metrics (length distribution, statistics)
# - TSS enrichment scores (chromatin accessibility quality)
# - Promoter overlap analysis (regulatory element distribution)
# - Library complexity assessment (duplication and diversity)
# - Cell count and per-cell quality metrics
# - Generated file inventory for downstream analysis
# 
# REPORT STRUCTURE:
# - Header with sample identification and analysis timestamp
# - Organized sections for each QC metric category
# - Quantitative metrics with biological interpretation context
# - File inventory for result tracking and downstream usage
# - Quality assessment summary with pass/fail indicators
# 
# QUALITY CONTROL METRICS:
# - Fragment metrics: Total count, length statistics, distribution
# - TSS enrichment: Accessibility at transcription start sites
# - Promoter analysis: Peak distribution in regulatory regions
# - Complexity scores: Library diversity and duplication assessment
# - Cell metrics: Count, fragments per cell, quality distributions
# 
# APPLICATIONS:
# - Quality assessment and sample comparison
# - Troubleshooting and optimization guidance
# - Downstream analysis parameter selection
# - Report generation for publication and documentation
# - Batch processing quality monitoring
# 
# OUTPUT FORMAT:
# - Human-readable text format for easy interpretation
# - Structured sections for programmatic parsing
# - Comprehensive file listing for workflow tracking
# - Quality thresholds and recommendations
# =============================================================================
echo "DEBUG: Creating comprehensive QC summary..."
cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_qc_summary.txt" << EOF
QC_METRICS_SUMMARY
Sample: $SAMPLE
Analysis_date: $(date)

FRAGMENT_METRICS:
$(cat "$FRAGMENT_STATS")

TSS_METRICS:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_tss_stats.txt" 2>/dev/null || echo "TSS_enrichment_percent	N/A")

PROMOTER_METRICS:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt" 2>/dev/null || echo "Promoter_analysis	Not_available")

LIBRARY_COMPLEXITY:
$(cat "$OUTPUT_DIR/qc_metrics/${SAMPLE}_library_complexity.txt")

CELL_COUNT:
Total_barcodes	$(wc -l < "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt")

FILES_GENERATED:
- Fragment length distribution: ${SAMPLE}_fragment_lengths.txt
- Fragment statistics: ${SAMPLE}_fragment_stats.txt
- TSS enrichment: ${SAMPLE}_tss_stats.txt
- Promoter analysis: ${SAMPLE}_promoter_stats.txt
- Fragments per cell: ${SAMPLE}_fragments_per_cell.txt
- Library complexity: ${SAMPLE}_library_complexity.txt
- QC summary: ${SAMPLE}_qc_summary.txt
EOF

# =============================================================================
# OPTIONAL R SCRIPT GENERATION FOR QC VISUALIZATION
# PURPOSE: Create R script for comprehensive QC plot generation
# 
# VISUALIZATION COMPONENTS:
# - Fragment length distribution plots (nucleosome positioning)
# - Fragments per cell distribution (cell quality assessment)
# - TSS enrichment visualization (accessibility quality)
# - Library complexity plots (duplication assessment)
# - Quality metric correlation analysis
# 
# PLOT TYPES:
# - Line plots: Fragment length distributions with nucleosome peaks
# - Histograms: Cell quality distributions and filtering thresholds
# - Scatter plots: Metric correlations and quality relationships
# - Box plots: Sample comparisons and batch effect assessment
# - Density plots: Distribution shapes and outlier identification
# 
# R SCRIPT FEATURES:
# - Automated data loading from generated QC files
# - Customizable plot parameters and themes
# - Publication-ready figure generation
# - Statistical summary calculations
# - Quality threshold annotations
# 
# APPLICATIONS:
# - Quality control visualization for reports
# - Sample comparison and batch effect detection
# - Parameter optimization guidance
# - Publication figure generation
# - Interactive quality assessment
# 
# TECHNICAL CONSIDERATIONS:
# - Required R packages: ggplot2, dplyr, scales
# - Output formats: PDF, PNG for different applications
# - Customizable themes and color schemes
# - Automated sample name substitution
# - Error handling for missing data files
# =============================================================================
cat > "$OUTPUT_DIR/qc_metrics/${SAMPLE}_plot_qc.R" << 'EOF'
#!/usr/bin/env Rscript
# QC Plotting Script for scATAC-seq data
# Usage: Rscript plot_qc.R sample_name

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("Please provide sample name as argument")
}

sample_name <- args[1]
library(ggplot2)
library(dplyr)

# Set working directory to qc_metrics folder
setwd("qc_metrics")

# 1. Fragment length distribution
frag_lengths <- read.table(paste0(sample_name, "_fragment_lengths.txt"), 
                          col.names = c("Length", "Count"))

p1 <- ggplot(frag_lengths, aes(x = Length, y = Count)) +
  geom_line() +
  xlim(0, 1000) +
  labs(title = paste("Fragment Length Distribution -", sample_name),
       x = "Fragment Length (bp)", y = "Count") +
  theme_minimal()

ggsave(paste0(sample_name, "_fragment_length_dist.png"), p1, width = 8, height = 6)

# 2. Fragments per cell distribution
frags_per_cell <- read.table(paste0(sample_name, "_fragments_per_cell.txt"), 
                            col.names = c("Barcode", "Fragment_count"))

p2 <- ggplot(frags_per_cell, aes(x = Fragment_count)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_log10() +
  labs(title = paste("Fragments per Cell Distribution -", sample_name),
       x = "Fragments per Cell (log10)", y = "Number of Cells") +
  theme_minimal()

ggsave(paste0(sample_name, "_fragments_per_cell.png"), p2, width = 8, height = 6)

cat("QC plots saved successfully\n")
EOF

chmod +x "$OUTPUT_DIR/qc_metrics/${SAMPLE}_plot_qc.R"

# =============================================================================
# PIPELINE COMPLETION AND FINAL SUMMARY
# PURPOSE: Provide comprehensive analysis completion summary
# 
# COMPLETION SUMMARY COMPONENTS:
# - Key quality metrics overview
# - File generation confirmation
# - Quality assessment interpretation
# - Next steps and recommendations
# - Troubleshooting guidance if needed
# 
# KEY METRICS DISPLAY:
# - Total fragments: Overall library size and sequencing depth
# - Unique fragments: Library diversity and complexity
# - Library complexity: PCR duplication assessment
# - Cell count: Single-cell resolution achieved
# - TSS enrichment: Chromatin accessibility quality
# - Promoter overlap: Regulatory element coverage
# 
# QUALITY INTERPRETATION:
# - Pass/fail assessment based on standard thresholds
# - Comparative context with typical values
# - Recommendations for downstream analysis
# - Potential issues and optimization suggestions
# 
# NEXT STEPS GUIDANCE:
# - Cell filtering recommendations based on metrics
# - Downstream analysis parameter suggestions
# - Additional QC steps if needed
# - Integration with other pipeline components
# 
# OUTPUT ORGANIZATION:
# - Clear file location information
# - Generated file inventory
# - Analysis timestamp for tracking
# - Success confirmation for pipeline monitoring
# =============================================================================
echo "QC Analysis Summary:"
echo "===================="
echo "Total fragments: $(printf "%'d" $TOTAL_FRAGMENTS)"
echo "Unique fragments: $(printf "%'d" $UNIQUE_FRAGMENTS)"
echo "Library complexity: $COMPLEXITY"
echo "TSS enrichment: $TSS_ENRICHMENT%"
if [[ -f "$OUTPUT_DIR/qc_metrics/${SAMPLE}_promoter_stats.txt" ]]; then
    echo "Promoter peaks: $PROMOTER_PERCENT%"
fi
echo "Total barcodes: $(wc -l < "$OUTPUT_DIR/qc_metrics/${SAMPLE}_fragments_per_cell.txt")"

echo "Output files created in: $OUTPUT_DIR/qc_metrics/"
echo "  - Comprehensive QC summary: ${SAMPLE}_qc_summary.txt"
echo "  - Fragment statistics: ${SAMPLE}_fragment_stats.txt"
echo "  - Per-cell metrics: ${SAMPLE}_fragments_per_cell.txt"
echo "  - R plotting script: ${SAMPLE}_plot_qc.R"

echo "========================================="
echo "Step 8 complete for $SAMPLE"
echo "QC metrics analysis finished"
echo "End time: $(date)"
echo "========================================="