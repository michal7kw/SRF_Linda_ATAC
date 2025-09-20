#!/bin/bash
# =============================================================================
# SCRIPT: 07_create_bigwig.sh
# PURPOSE: Generate BigWig files for ATAC-seq data visualization
# 
# DESCRIPTION:
# This script converts aligned ATAC-seq reads into BigWig format files for
# genome browser visualization and downstream analysis. Creates both raw
# coverage and RPM-normalized tracks for comprehensive data exploration.
# 
# INPUT REQUIREMENTS:
# - Aligned reads in BED format (from chromap alignment)
# - Chromosome sizes file for genome reference
# - Properly configured conda environment with visualization tools
# 
# OUTPUT FILES:
# - Raw coverage BigWig: Direct read coverage tracks
# - RPM-normalized BigWig: Reads per million normalized tracks
# - Summary statistics: Processing metadata and file information
# 
# COMPUTATIONAL REQUIREMENTS:
# - CPU: 4 cores for parallel processing
# - Memory: 16GB for large genome coverage calculations
# - Time: Up to 4 hours for comprehensive datasets
# - Storage: Significant space for intermediate and final files
# 
# DEPENDENCIES:
# - bedtools: Genome coverage calculation
# - UCSC tools: bedGraphToBigWig conversion
# - bc calculator: Normalization factor computation
# =============================================================================
#SBATCH --job-name=create_bigwig
#SBATCH --output=logs/07_create_bigwig_%a.out
#SBATCH --error=logs/07_create_bigwig_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=workq

# =============================================================================
# CONDA ENVIRONMENT SETUP
# PURPOSE: Initialize specialized environment for BigWig file generation
# 
# ENVIRONMENT FEATURES:
# - bedtools: Genome coverage and interval operations
# - UCSC tools: bedGraphToBigWig conversion utilities
# - bc calculator: Mathematical operations for normalization
# - Standard Unix tools: File processing and validation
# 
# ERROR HANDLING:
# - Strict mode enabled (set -euo pipefail)
# - Immediate exit on command failures
# - Undefined variable detection
# - Pipeline failure propagation
# 
# VISUALIZATION REQUIREMENTS:
# - High-performance tools for large-scale data processing
# - Memory-efficient algorithms for genome-wide coverage
# - Standardized output formats for browser compatibility
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

set -euo pipefail

# =============================================================================
# CONFIGURATION AND SAMPLE SETUP
# PURPOSE: Define sample parameters and processing environment
# 
# SAMPLE CONFIGURATION:
# - Array-based processing for parallel execution
# - Standardized naming convention for consistency
# - Automated sample selection via SLURM array indexing
# 
# PATH CONFIGURATION:
# - Centralized output directory for organized results
# - Consistent file structure across pipeline steps
# - Shared access for downstream analysis tools
# 
# ARRAY PROCESSING:
# - Enables parallel BigWig generation across samples
# - Optimizes cluster resource utilization
# - Maintains independent processing for each sample
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 7: Creating BigWig files for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY PREPARATION
# PURPOSE: Ensure organized directory structure for BigWig files
# 
# DIRECTORY STRUCTURE:
# - bigwig/: Primary BigWig files and metadata
# - Separate namespace for visualization files
# - Organized storage for browser track files
# 
# BENEFITS:
# - Clear separation from analysis results
# - Easy identification of visualization files
# - Simplified file management and sharing
# =============================================================================
mkdir -p "$OUTPUT_DIR/bigwig"

# =============================================================================
# PREREQUISITE VALIDATION
# PURPOSE: Verify required input files from previous pipeline steps
# 
# VALIDATION CHECKS:
# - Reads file existence from chromap alignment
# - Chromosome sizes file for genome reference
# - File accessibility and readability
# 
# INPUT REQUIREMENTS:
# - BED format reads with proper coordinates
# - Chromosome sizes in standard UCSC format
# - Compressed files for storage efficiency
# 
# ERROR HANDLING:
# - Clear error messages for missing dependencies
# - Guidance for resolving prerequisite issues
# - Graceful exit with informative feedback
# =============================================================================
READS_FILE="$OUTPUT_DIR/${SAMPLE}_reads.bed.gz"
CHROM_SIZES="$OUTPUT_DIR/mm10.chrom.sizes"

if [[ ! -f "$READS_FILE" ]]; then
    echo "ERROR: Reads file not found: $READS_FILE"
    echo "Please run chromap alignment first"
    exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "ERROR: Chromosome sizes file not found: $CHROM_SIZES"
    echo "Please ensure mm10.chrom.sizes exists in output directory"
    exit 1
fi

# =============================================================================
# TOOL VALIDATION AND VERIFICATION
# PURPOSE: Ensure all required tools are available for BigWig generation
# 
# REQUIRED TOOLS:
# - bedtools: Genome coverage calculation from BED files
# - bedGraphToBigWig: UCSC tool for bedGraph to BigWig conversion
# - bc calculator: Mathematical operations for normalization
# 
# TOOL FUNCTIONS:
# - bedtools genomecov: Creates coverage tracks from read positions
# - bedGraphToBigWig: Converts text bedGraph to binary BigWig format
# - bc: Calculates RPM normalization factors
# 
# INSTALLATION GUIDANCE:
# - Provides specific conda installation commands
# - Ensures reproducible environment setup
# - Validates tool accessibility and functionality
# 
# ERROR HANDLING:
# - Clear error messages for missing tools
# - Installation instructions for resolution
# - Graceful exit with informative feedback
# =============================================================================
echo "DEBUG: Checking for required tools..."
if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools not found. Please install bedtools"
    echo "DEBUG: You can install with: conda install -c bioconda bedtools"
    exit 1
fi

if ! command -v bedGraphToBigWig &> /dev/null; then
    echo "ERROR: bedGraphToBigWig not found. Please install UCSC tools"
    echo "DEBUG: You can install with: conda install -c bioconda ucsc-bedgraphtobigwig"
    exit 1
fi

echo "DEBUG: bedtools found: $(which bedtools)"
echo "DEBUG: bedGraphToBigWig found: $(which bedGraphToBigWig)"

# =============================================================================
# RAW COVERAGE BIGWIG GENERATION
# PURPOSE: Create genome-wide coverage tracks from aligned reads
# 
# PROCESSING WORKFLOW:
# 1. Generate bedGraph from BED reads using bedtools genomecov
# 2. Convert bedGraph to BigWig using UCSC tools
# 3. Validate output and clean up intermediate files
# 
# BEDGRAPH FORMAT:
# chr1    1000    1200    5
# chr1    1200    1400    3
# (chromosome, start, end, coverage_value)
# 
# BIGWIG ADVANTAGES:
# - Binary format for efficient storage and access
# - Indexed for fast random access
# - Compatible with genome browsers (IGV, UCSC)
# - Supports zooming and panning operations
# 
# COVERAGE CALCULATION:
# - Uses bedtools genomecov with -bg flag
# - Generates coverage for all genomic positions
# - Excludes zero-coverage regions for efficiency
# - Maintains strand-independent coverage
# 
# FILE MANAGEMENT:
# - Creates intermediate bedGraph for conversion
# - Removes bedGraph after successful BigWig creation
# - Validates file sizes and content
# - Skips processing if BigWig already exists
# =============================================================================
BEDGRAPH_FILE="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bedGraph"
BIGWIG_FILE="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bw"

if [[ ! -f "$BIGWIG_FILE" ]]; then
    echo "DEBUG: Creating coverage bedGraph..."
    echo "DEBUG: This may take several minutes for large datasets..."
    
    # Create bedGraph with coverage information
    # Using bedtools genomecov to create coverage track from BED file
    zcat "$READS_FILE" | \
        bedtools genomecov -i stdin -g "$CHROM_SIZES" -bg > "$BEDGRAPH_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: BedGraph created successfully"
        echo "DEBUG: BedGraph file size: $(stat -c%s "$BEDGRAPH_FILE") bytes"
        echo "DEBUG: Number of intervals: $(wc -l < "$BEDGRAPH_FILE")"
    else
        echo "ERROR: Failed to create bedGraph"
        exit 1
    fi
    
    # Convert bedGraph to BigWig
    echo "DEBUG: Converting bedGraph to BigWig..."
    bedGraphToBigWig "$BEDGRAPH_FILE" "$CHROM_SIZES" "$BIGWIG_FILE"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: BigWig created successfully"
        echo "DEBUG: BigWig file size: $(stat -c%s "$BIGWIG_FILE") bytes"
        
        # Clean up intermediate bedGraph file
        rm -f "$BEDGRAPH_FILE"
        echo "DEBUG: Cleaned up intermediate bedGraph file"
    else
        echo "ERROR: Failed to convert bedGraph to BigWig"
        exit 1
    fi
else
    echo "DEBUG: Skipping BigWig creation, file already exists: $BIGWIG_FILE"
fi

# =============================================================================
# RPM-NORMALIZED BIGWIG GENERATION
# PURPOSE: Create reads-per-million normalized coverage tracks
# 
# RPM NORMALIZATION:
# - Normalizes coverage by total read count
# - Formula: (coverage * 1,000,000) / total_reads
# - Enables comparison between samples with different sequencing depths
# - Standard normalization for ATAC-seq visualization
# 
# NORMALIZATION WORKFLOW:
# 1. Count total reads in the dataset
# 2. Calculate RPM scaling factor (1M / total_reads)
# 3. Apply scaling during bedGraph generation
# 4. Convert normalized bedGraph to BigWig
# 
# BIOLOGICAL SIGNIFICANCE:
# - Accounts for sequencing depth differences
# - Enables direct comparison between samples
# - Reveals relative accessibility patterns
# - Standard for publication-quality figures
# 
# COMPUTATIONAL CONSIDERATIONS:
# - Uses bc calculator for precise floating-point arithmetic
# - Applies scaling factor during bedtools genomecov
# - Maintains precision with 6 decimal places
# - Efficient single-pass processing
# 
# OUTPUT APPLICATIONS:
# - Comparative analysis between samples
# - Meta-analysis and data integration
# - Publication figures and heatmaps
# - Quantitative accessibility measurements
# =============================================================================
NORMALIZED_BIGWIG="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bw"
NORMALIZED_BEDGRAPH="$OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bedGraph"

if [[ ! -f "$NORMALIZED_BIGWIG" ]]; then
    echo "DEBUG: Creating RPM-normalized BigWig..."
    
    # Count total reads for normalization
    TOTAL_READS=$(zcat "$READS_FILE" | wc -l)
    SCALE_FACTOR=$(echo "scale=6; 1000000 / $TOTAL_READS" | bc)
    
    echo "DEBUG: Total reads: $(printf "%'d" $TOTAL_READS)"
    echo "DEBUG: RPM scale factor: $SCALE_FACTOR"
    
    # Create RPM-normalized bedGraph
    zcat "$READS_FILE" | \
        bedtools genomecov -i stdin -g "$CHROM_SIZES" -bg -scale $SCALE_FACTOR > "$NORMALIZED_BEDGRAPH"
    
    if [[ $? -eq 0 ]]; then
        echo "DEBUG: RPM-normalized bedGraph created"
        
        # Convert to BigWig
        bedGraphToBigWig "$NORMALIZED_BEDGRAPH" "$CHROM_SIZES" "$NORMALIZED_BIGWIG"
        
        if [[ $? -eq 0 ]]; then
            echo "DEBUG: RPM-normalized BigWig created successfully"
            rm -f "$NORMALIZED_BEDGRAPH"
            echo "DEBUG: Cleaned up intermediate RPM bedGraph file"
        else
            echo "ERROR: Failed to convert RPM bedGraph to BigWig"
            exit 1
        fi
    else
        echo "ERROR: Failed to create RPM-normalized bedGraph"
        exit 1
    fi
else
    echo "DEBUG: Skipping RPM-normalized BigWig creation, file already exists: $NORMALIZED_BIGWIG"
fi

# =============================================================================
# SUMMARY STATISTICS AND METADATA GENERATION
# PURPOSE: Document BigWig generation process and output characteristics
# 
# SUMMARY COMPONENTS:
# - Processing metadata (sample name, timestamp)
# - Input file documentation (reads, chromosome sizes)
# - Output file inventory (raw and normalized BigWig)
# - Quantitative statistics (read counts, file sizes)
# - Processing status and validation
# 
# METADATA IMPORTANCE:
# - Enables reproducibility and traceability
# - Documents processing parameters and inputs
# - Facilitates quality control and validation
# - Supports downstream analysis planning
# 
# FILE SIZE ANALYSIS:
# - Indicates compression efficiency
# - Helps estimate storage requirements
# - Validates successful file creation
# - Enables performance optimization
# 
# QUALITY CONTROL METRICS:
# - Total read count validation
# - File existence verification
# - Size consistency checks
# - Processing completion confirmation
# 
# DOWNSTREAM APPLICATIONS:
# - Analysis pipeline documentation
# - Batch processing summaries
# - Quality control reports
# - Method documentation for publications
# =============================================================================
echo "DEBUG: Creating BigWig summary statistics..."
cat > "$OUTPUT_DIR/bigwig/${SAMPLE}_bigwig_info.txt" << EOF
SAMPLE=$SAMPLE
BIGWIG_CREATED_AT=$(date)
INPUT_READS_FILE=$READS_FILE
TOTAL_READS=$TOTAL_READS
RPM_SCALE_FACTOR=$SCALE_FACTOR
RAW_BIGWIG_FILE=${SAMPLE}_coverage.bw
NORMALIZED_BIGWIG_FILE=${SAMPLE}_coverage_rpm.bw
RAW_BIGWIG_SIZE=$(stat -c%s "$BIGWIG_FILE")
NORMALIZED_BIGWIG_SIZE=$(stat -c%s "$NORMALIZED_BIGWIG")
EOF

# =============================================================================
# PIPELINE COMPLETION AND FINAL SUMMARY
# PURPOSE: Confirm successful BigWig generation and provide final status
# 
# COMPLETION VERIFICATION:
# - Confirms successful processing of sample
# - Reports total reads processed
# - Documents output location
# - Provides completion timestamp
# 
# SUCCESS INDICATORS:
# - Sample name confirmation
# - Read count validation
# - Output directory specification
# - Timestamp for processing duration
# 
# NEXT STEPS:
# - BigWig files ready for genome browser visualization
# - Coverage tracks available for downstream analysis
# - Files suitable for comparative studies
# - Ready for integration with other pipeline outputs
# 
# QUALITY ASSURANCE:
# - Validates complete pipeline execution
# - Confirms output file generation
# - Documents processing completion
# - Enables batch processing verification
# =============================================================================
echo "BigWig files created:"
echo "  - Raw coverage: $OUTPUT_DIR/bigwig/${SAMPLE}_coverage.bw"
echo "  - RPM normalized: $OUTPUT_DIR/bigwig/${SAMPLE}_coverage_rpm.bw"
echo "  - Summary info: $OUTPUT_DIR/bigwig/${SAMPLE}_bigwig_info.txt"

echo "========================================="
echo "Step 7 complete for $SAMPLE"
echo "Created BigWig files for $(printf "%'d" $TOTAL_READS) reads"
echo "End time: $(date)"
echo "========================================="