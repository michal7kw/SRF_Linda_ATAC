#!/bin/bash
# =============================================================================
# SCRIPT: 04_chromap_alignment.sh
# PURPOSE: Perform chromap alignment for scATAC-seq data with barcode integration
# 
# DESCRIPTION:
# This script uses chromap to align paired-end genomic reads (R1/R3) while
# incorporating cell barcode information from extracted R2 reads. Chromap is
# optimized for single-cell ATAC-seq data, providing efficient alignment and
# fragment generation with barcode-aware processing.
# 
# KEY OPERATIONS:
# 1. Load optimal barcode parameters from step 2 testing results
# 2. Install/configure chromap alignment tool
# 3. Build or verify chromap genome index
# 4. Perform barcode-aware alignment with error correction
# 5. Generate compressed and indexed fragment files
# 6. Calculate alignment statistics and quality metrics
# 
# INPUT REQUIREMENTS:
# - Validated R1, R3 genomic FASTQ files
# - Extracted and validated barcode file from steps 1-3
# - Barcode test results with optimal whitelist selection
# - Reference genome FASTA file
# 
# OUTPUT:
# - Fragment file (BED-like format with barcodes)
# - Alignment statistics and quality metrics
# - Compressed and indexed fragment files for downstream analysis
# 
# COMPUTATIONAL REQUIREMENTS:
# - High memory (128GB) for genome index and alignment
# - Multi-threading (16 cores) for efficient processing
# - Extended runtime (24h) for large datasets
# =============================================================================
#SBATCH --job-name=chromap_align
#SBATCH --output=logs/04_chromap_align_%a.out
#SBATCH --error=logs/04_chromap_align_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=workq

# =============================================================================
# ENVIRONMENT SETUP
# PURPOSE: Initialize conda environment and error handling for alignment
# 
# OPERATIONS:
# - Activates alignment_two conda environment
# - Enables strict error handling for robust pipeline execution
# 
# REQUIRED TOOLS:
# - chromap: Single-cell ATAC-seq alignment tool
# - bgzip/tabix: Fragment file compression and indexing
# - Standard Unix utilities for file processing
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# ALIGNMENT CONFIGURATION
# PURPOSE: Define samples, paths, and alignment parameters
# 
# SAMPLE PROCESSING:
# - Uses SLURM array for parallel sample processing
# - Each job processes one sample independently
# 
# DIRECTORY STRUCTURE:
# - DATA_DIR: Original FASTQ files
# - PROCESSED_DATA_DIR: Extracted barcode files
# - OUTPUT_DIR: Alignment results and fragment files
# 
# REFERENCE GENOME:
# - mm10.fa: Mouse reference genome for alignment
# - CHROMAP_INDEX: Pre-built or generated alignment index
# 
# PERFORMANCE PARAMETERS:
# - THREADS=16: Parallel processing for alignment efficiency
# - Optimized for high-memory, multi-core environment
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/chromap_final_output"
REF_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"
CHROMAP_INDEX="$OUTPUT_DIR/chromap_index/mm10.index"
THREADS=16

echo "========================================="
echo "Step 4: Chromap alignment for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY PREPARATION
# PURPOSE: Create necessary directories for alignment outputs
# 
# DIRECTORIES:
# - fragments/: Fragment files and indices
# - logs/: Alignment logs and statistics
# =============================================================================
mkdir -p "$OUTPUT_DIR/fragments" "$OUTPUT_DIR/logs"

EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

# =============================================================================
# PREREQUISITE VALIDATION
# PURPOSE: Verify that all required input files exist before alignment
# 
# DEPENDENCIES:
# - Extracted barcode file from steps 1-3 (count validation)
# - Barcode test results from step 2 (optimal whitelist selection)
# 
# ERROR HANDLING:
# - Clear error messages with guidance on missing prerequisites
# - Prevents alignment with incomplete or invalid inputs
# =============================================================================
if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run steps 1-3 first"
    exit 1
fi

# Load barcode test results
BARCODE_RESULTS="$OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt"
if [[ ! -f "$BARCODE_RESULTS" ]]; then
    echo "ERROR: Barcode test results not found: $BARCODE_RESULTS"
    echo "Please run step 2 (02_test_barcodes.sh) first"
    exit 1
fi

# =============================================================================
# BARCODE CONFIGURATION LOADING
# PURPOSE: Load optimal barcode parameters determined in step 2
# 
# LOADED VARIABLES:
# - BEST_WHITELIST: Optimal whitelist orientation (forward/reverse)
# - MATCH_RATE: Percentage of barcodes matching the selected whitelist
# - BC_ERROR_THRESHOLD: Error tolerance for barcode matching
# 
# RATIONALE:
# Step 2 testing determines the best whitelist orientation and error threshold
# based on actual data characteristics, ensuring optimal barcode matching
# during alignment.
# =============================================================================
source "$BARCODE_RESULTS"

echo "DEBUG: Using barcode test results:"
echo "  Best whitelist: $BEST_WHITELIST"
echo "  Match rate: $MATCH_RATE%"
echo "  Error threshold: $BC_ERROR_THRESHOLD"

# =============================================================================
# CHROMAP INSTALLATION AND VERIFICATION
# PURPOSE: Ensure chromap alignment tool is available and functional
# 
# INSTALLATION STRATEGY:
# 1. Check if chromap is already available in PATH
# 2. Attempt to load system module with multiple naming conventions
# 3. Fall back to conda installation if module unavailable
# 4. Verify successful installation before proceeding
# 
# ERROR HANDLING:
# - Graceful fallback between installation methods
# - Clear error messages with manual installation guidance
# - Version verification for debugging and reproducibility
# 
# SUPPORTED ENVIRONMENTS:
# - HPC clusters with module systems
# - Conda-based environments
# - Manual installations in PATH
# =============================================================================
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

# =============================================================================
# CHROMAP INDEX PREPARATION
# PURPOSE: Create or verify chromap genome index for efficient alignment
# 
# INDEX FUNCTIONALITY:
# - Pre-computed data structure for fast genome sequence lookup
# - Enables rapid alignment of millions of reads
# - Shared across samples to avoid redundant computation
# 
# OPERATIONS:
# 1. Check if index already exists (reuse for efficiency)
# 2. Create index directory if needed
# 3. Build index from reference genome FASTA
# 4. Verify successful index creation
# 
# PERFORMANCE IMPACT:
# - Index creation: One-time cost (~30-60 minutes for mouse genome)
# - Index reuse: Significant time savings for subsequent runs
# - Memory requirement: ~8-16GB during index creation
# - Disk space: ~4-8GB for completed index
# =============================================================================
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

# =============================================================================
# CHROMAP ALIGNMENT EXECUTION
# PURPOSE: Perform barcode-aware alignment of scATAC-seq reads
# 
# ALIGNMENT PARAMETERS:
# --preset atac: Optimized settings for ATAC-seq data
# -x: Pre-built genome index for fast alignment
# -r: Reference genome FASTA file
# -1/-2: Paired-end genomic reads (R1/R3)
# -b: Extracted barcode file (R2)
# --barcode-whitelist: Validated barcode whitelist from step 2
# --bc-error-threshold: Error tolerance for barcode matching
# -e: Maximum edit distance for alignment
# -t: Number of threads for parallel processing
# --low-mem: Memory-efficient mode for large datasets
# -o: Output fragment file path
# 
# OUTPUT FORMAT:
# Fragment file contains: chr, start, end, barcode, read_count
# =============================================================================
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" ]]; then
    echo "DEBUG: Running chromap with the following parameters:"
    echo "  Reference: $REF_GENOME"
    echo "  Index: $CHROMAP_INDEX"
    echo "  R1: $DATA_DIR/${SAMPLE}_R1_001.fastq.gz"
    echo "  R2: $DATA_DIR/${SAMPLE}_R3_001.fastq.gz"
    echo "  Barcode: $EXTRACTED_BC_FILE (rightmost 16bp from R2)"
    echo "  Whitelist: $WHITELIST"
    echo "  Barcode error threshold: $BC_ERROR_THRESHOLD"
    echo "  Expected match rate: ${MATCH_RATE}%"

    chromap --preset atac \
        -x "$CHROMAP_INDEX" \
        -r "$REF_GENOME" \
        -1 "$DATA_DIR/${SAMPLE}_R1_001.fastq.gz" \
        -2 "$DATA_DIR/${SAMPLE}_R3_001.fastq.gz" \
        -b "$EXTRACTED_BC_FILE" \
        --barcode-whitelist "$WHITELIST" \
        --bc-error-threshold 2 \
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

# =============================================================================
# FRAGMENT FILE COMPRESSION AND INDEXING
# PURPOSE: Optimize fragment files for efficient downstream analysis
# 
# COMPRESSION BENEFITS:
# - Reduces file size by ~70-80% (important for large datasets)
# - Maintains random access capability through indexing
# - Standard format for scATAC-seq analysis tools
# 
# INDEXING FUNCTIONALITY:
# - Creates tabix index for rapid genomic region queries
# - Enables efficient overlap detection with genomic features
# - Required for most downstream analysis tools
# 
# TOOLS USED:
# - bgzip: Block-compressed gzip (maintains random access)
# - tabix: Creates coordinate-sorted index for genomic intervals
# =============================================================================
echo "DEBUG: Processing fragments..."
if [[ ! -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz.tbi" ]]; then
    bgzip -f "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv"
    tabix -f -s 1 -b 2 -e 3 -p bed "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz"
    echo "DEBUG: Fragments compressed and indexed"
else
    echo "DEBUG: Skipping fragment indexing, index file already exists."
fi

# =============================================================================
# FRAGMENT STATISTICS GENERATION
# PURPOSE: Calculate comprehensive alignment metrics for quality assessment
# 
# METRICS CALCULATED:
# - Total fragment count: Overall alignment success indicator
# - Unique barcodes: Cell diversity and barcode quality
# - Fragments in cells: Fraction matching validated whitelist
# - Fraction in cells: Quality metric for barcode assignment
# 
# USAGE:
# - Quality control for alignment success
# - Comparison between samples
# - Input for downstream analysis planning
# =============================================================================
TOTAL_FRAGS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | wc -l)
UNIQUE_BC=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | sort -u | wc -l)
FRAGS_IN_CELLS=$(zcat "$OUTPUT_DIR/fragments/${SAMPLE}_fragments.tsv.gz" | cut -f4 | grep -Ff "$WHITELIST" | wc -l)

echo "Fragment statistics:"
echo "  Total fragments: $TOTAL_FRAGS"
echo "  Unique barcodes: $UNIQUE_BC"
echo "  Fragments in cells: $FRAGS_IN_CELLS"
echo "  Fraction in cells: $(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)"

# =============================================================================
# ALIGNMENT RESULTS PERSISTENCE
# PURPOSE: Save comprehensive alignment summary for pipeline tracking
# 
# RECORDED INFORMATION:
# - Sample identification and processing timestamp
# - Fragment statistics and quality metrics
# - Barcode configuration used (whitelist and error threshold)
# - Alignment success confirmation
# 
# USAGE:
# - Pipeline progress tracking
# - Parameter documentation for reproducibility
# - Quality control and troubleshooting reference
# =============================================================================
cat > "$OUTPUT_DIR/logs/${SAMPLE}_alignment_results.txt" << EOF
SAMPLE=$SAMPLE
TOTAL_FRAGS=$TOTAL_FRAGS
UNIQUE_BC=$UNIQUE_BC
FRAGS_IN_CELLS=$FRAGS_IN_CELLS
FRACTION_IN_CELLS=$(echo "scale=4; $FRAGS_IN_CELLS / $TOTAL_FRAGS" | bc)
WHITELIST_USED=$WHITELIST
MATCH_RATE_USED=$MATCH_RATE
BC_ERROR_THRESHOLD_USED=$BC_ERROR_THRESHOLD
ALIGNMENT_COMPLETED_AT=$(date)
EOF

echo "Results saved to: $OUTPUT_DIR/logs/${SAMPLE}_alignment_results.txt"

echo "========================================="
echo "Step 4 complete for $SAMPLE"
echo "Generated $(printf "%'d" $TOTAL_FRAGS) fragments from $(printf "%'d" $UNIQUE_BC) unique barcodes"
echo "End time: $(date)"
echo "========================================="