#!/bin/bash
# =============================================================================
# SCRIPT: 02_test_barcodes.sh
# PURPOSE: Test extracted barcodes against 10x Genomics whitelists for quality validation
# 
# FUNCTIONALITY:
# - Downloads and prepares 10x Genomics ATAC barcode whitelists
# - Creates reverse complement whitelist for orientation testing
# - Tests extracted barcodes against both whitelist orientations
# - Determines optimal whitelist and barcode error correction threshold
# - Validates barcode length and N-content
# 
# INPUT: Extracted 16bp barcode FASTQ files from step 1
# OUTPUT: Barcode test results and optimal whitelist selection
# 
# EDGE CASES:
# - Handles whitelist download failures
# - Implements timeout protection for large file operations
# - Provides fallback matching algorithms for performance
# - Adjusts error thresholds based on match rates
# =============================================================================

#SBATCH --job-name=test_barcodes
#SBATCH --output=logs/02_test_barcodes_%a.out
#SBATCH --error=logs/02_test_barcodes_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# =============================================================================
# ENVIRONMENT SETUP
# PURPOSE: Initialize conda environment and set strict error handling
# 
# OPERATIONS:
# - Activates conda environment with required bioinformatics tools
# - Enables strict bash error handling for robust execution
# =============================================================================
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate alignment_two

set -euo pipefail

# =============================================================================
# CONFIGURATION AND DIRECTORY SETUP
# PURPOSE: Define sample names and directory paths for barcode testing
# 
# OPERATIONS:
# - Maps SLURM array task ID to specific sample name
# - Sets input directory containing extracted barcode files
# - Sets output directory for QC results and whitelist files
# 
# INPUT EXPECTATIONS:
# - Extracted barcode files from step 1 in PROCESSED_DATA_DIR
# - SLURM_ARRAY_TASK_ID: 0 or 1 (corresponding to control/mutant samples)
# =============================================================================
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

PROCESSED_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/ATAC_data/DATA/nestin_processed"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/ATAC_data/chromap_final_output"

echo "========================================="
echo "Step 2: Testing barcode quality for $SAMPLE"
echo "Start time: $(date)"
echo "========================================="

# =============================================================================
# OUTPUT DIRECTORY CREATION AND FILE VALIDATION
# PURPOSE: Ensure required directories exist and validate input files
# 
# OPERATIONS:
# - Creates QC output directory for storing test results
# - Validates that extracted barcode file from step 1 exists
# 
# ERROR HANDLING:
# - Exits with clear error message if prerequisite files are missing
# - Prevents execution without proper input files
# =============================================================================
mkdir -p "$OUTPUT_DIR/qc"

EXTRACTED_BC_FILE="$PROCESSED_DATA_DIR/${SAMPLE}_R2_rightmost16bp.fastq.gz"

if [[ ! -f "$EXTRACTED_BC_FILE" ]]; then
    echo "ERROR: Extracted barcode file not found: $EXTRACTED_BC_FILE"
    echo "Please run step 1 (01_extract_barcodes.sh) first"
    exit 1
fi

# =============================================================================
# WHITELIST DOWNLOAD AND PREPARATION
# PURPOSE: Download and prepare 10x Genomics ATAC barcode whitelists
# 
# OPERATIONS:
# - Creates whitelist directory if it doesn't exist
# - Downloads official 10x Genomics ATAC 737K barcode whitelist
# - Extracts compressed whitelist file
# 
# WHITELIST DETAILS:
# - Source: 10x Genomics official ATAC barcode whitelist (737K barcodes)
# - Format: One barcode per line, 16bp sequences
# - Used for validating extracted barcodes against known valid sequences
# 
# ERROR HANDLING:
# - Validates successful download before proceeding
# - Exits with error if whitelist cannot be obtained
# =============================================================================
WHITELIST_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/ATAC_data/barcode_whitelists"
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

# =============================================================================
# REVERSE COMPLEMENT WHITELIST GENERATION
# PURPOSE: Create reverse complement whitelist for orientation testing
# 
# OPERATIONS:
# - Counts barcodes in original whitelist for validation
# - Generates reverse complement of each barcode using rev and tr commands
# - Creates separate RC whitelist file for comparison testing
# 
# ALGORITHM:
# - rev: Reverses each barcode sequence
# - tr 'ACGT' 'TGCA': Complements nucleotides (A↔T, C↔G)
# - Combined effect: Creates reverse complement sequences
# 
# RATIONALE:
# - Some sequencing protocols may produce barcodes in reverse complement
# - Testing both orientations ensures optimal barcode matching
# - Determines which orientation provides better match rates
# =============================================================================
echo "DEBUG: Checking whitelist files..."
echo "  Original whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1.txt") barcodes"

if [[ ! -f "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt" ]]; then
    echo "DEBUG: Creating reverse complement whitelist..."
    cat "$WHITELIST_DIR/atac_737K-arc-v1.txt" | rev | tr 'ACGT' 'TGCA' > \
        "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
fi

echo "  RC whitelist: $(wc -l < "$WHITELIST_DIR/atac_737K-arc-v1_rc.txt") barcodes"

# =============================================================================
# TEST BARCODE EXTRACTION AND VALIDATION
# PURPOSE: Extract subset of barcodes for whitelist testing and validation
# 
# OPERATIONS:
# - Extracts sequence lines (every 4th line starting from line 2) from FASTQ
# - Limits to first 10,000 barcodes for efficient testing
# - Implements timeout protection to prevent hanging on large files
# - Validates successful extraction and non-zero barcode count
# 
# ALGORITHM:
# - awk 'NR%4==2': Extracts only sequence lines from FASTQ format
# - head -10000: Limits to manageable subset for testing
# - timeout 300: Prevents indefinite hanging (5-minute limit)
# 
# ERROR HANDLING:
# - Detects timeout conditions (exit code 124)
# - Validates extraction success and non-empty output
# - Prevents downstream analysis with invalid data
# =============================================================================
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

# =============================================================================
# BARCODE LENGTH AND QUALITY VALIDATION
# PURPOSE: Validate barcode format and assess sequence quality
# 
# OPERATIONS:
# - Extracts sample barcode for length verification
# - Validates that barcodes are exactly 16bp as expected
# - Analyzes N-content to assess sequence quality
# 
# VALIDATION CHECKS:
# - Length validation: Critical for downstream compatibility
# - N-content analysis: High N-content indicates poor sequencing quality
# 
# ERROR HANDLING:
# - Exits if barcode length is incorrect (indicates extraction problems)
# - Reports N-content for quality assessment
# =============================================================================
SAMPLE_SEQ=$(head -1 "$OUTPUT_DIR/qc/test_barcodes.txt")
SEQ_LENGTH=${#SAMPLE_SEQ}
echo "DEBUG: Sample barcode: $SAMPLE_SEQ"
echo "DEBUG: Barcode length: $SEQ_LENGTH"

if [[ $SEQ_LENGTH -ne 16 ]]; then
    echo "ERROR: Expected 16bp barcodes, got ${SEQ_LENGTH}bp"
    exit 1
fi

# Analyze N content
N_COUNT=$(echo "$SAMPLE_SEQ" | grep -o 'N' | wc -l)
echo "DEBUG: Number of N's in sample barcode: $N_COUNT"

# =============================================================================
# WHITELIST ORIENTATION TESTING
# PURPOSE: Test both whitelist orientations to determine optimal matching
# 
# OPERATIONS:
# - Creates smaller test subset (1000 barcodes) for efficient testing
# - Tests barcodes against original whitelist using grep
# - Implements timeout protection and fallback algorithms
# - Uses both grep and awk methods for robustness
# 
# ALGORITHM DETAILS:
# - grep -Fxf: Fixed string exact match from file (fast but memory-intensive)
# - awk method: Hash-based lookup (slower but more memory-efficient)
# - timeout 60/30: Progressive timeout reduction for fallback methods
# 
# PERFORMANCE OPTIMIZATION:
# - Uses smaller subset (1000 vs 10000) for orientation testing
# - Implements multiple matching strategies to handle large whitelists
# - Provides fallback when primary method times out
# =============================================================================
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

# =============================================================================
# REVERSE COMPLEMENT WHITELIST TESTING
# PURPOSE: Test barcodes against reverse complement whitelist for comparison
# 
# OPERATIONS:
# - Tests same barcode subset against reverse complement whitelist
# - Uses identical timeout and fallback strategy as original testing
# - Calculates match percentages for both orientations
# 
# COMPARISON METRICS:
# - Absolute matches: Raw count of matching barcodes
# - Match percentage: Proportion of test barcodes matching whitelist
# - Used to determine which orientation provides better coverage
# =============================================================================
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

# =============================================================================
# OPTIMAL WHITELIST SELECTION
# PURPOSE: Select whitelist orientation with highest match rate
# 
# DECISION LOGIC:
# - Compares absolute match counts between orientations
# - Selects orientation with more matching barcodes
# - Sets variables for downstream processing
# 
# OUTPUT:
# - WHITELIST: Path to selected whitelist file
# - BEST_WHITELIST: Orientation identifier for reporting
# =============================================================================
if [[ $MATCHES_RC -gt $MATCHES_ORIG ]]; then
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1_rc.txt"
    BEST_WHITELIST="reverse_complement"
    echo "Using reverse complement whitelist"
else
    WHITELIST="$WHITELIST_DIR/atac_737K-arc-v1.txt"
    BEST_WHITELIST="original"
    echo "Using original whitelist"
fi

# =============================================================================
# BARCODE ERROR THRESHOLD DETERMINATION
# PURPOSE: Set error correction threshold based on observed match rates
# 
# ALGORITHM:
# - Calculates best match rate from both orientations
# - Sets error correction threshold based on match quality
# - Higher thresholds for lower match rates (more aggressive correction)
# 
# THRESHOLD LOGIC:
# - <5%: Critical error - likely fundamental problem
# - 5-20%: Low quality - aggressive error correction (threshold=2)
# - 20-50%: Moderate quality - standard error correction (threshold=1)
# - >50%: Good quality - minimal error correction (threshold=1)
# 
# ERROR HANDLING:
# - Exits if match rate is critically low (<5%)
# - Provides warnings for suboptimal match rates
# - Adjusts downstream processing parameters accordingly
# =============================================================================
MAX_MATCHES=$((MATCHES_ORIG > MATCHES_RC ? MATCHES_ORIG : MATCHES_RC))
MATCH_RATE=$(echo "scale=2; ${MAX_MATCHES:-0} * 100 / $SMALL_BC_COUNT" | bc)

echo "DEBUG: Best match rate: ${MATCH_RATE}% (from $SMALL_BC_COUNT test barcodes)"

if [[ $(echo "$MATCH_RATE < 5.0" | bc) -eq 1 ]]; then
    echo "ERROR: Very low barcode match rate (${MATCH_RATE}%)"
    BC_ERROR_THRESHOLD=3
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

# =============================================================================
# RESULTS PERSISTENCE
# PURPOSE: Save barcode testing results for downstream processing
# 
# OPERATIONS:
# - Creates structured results file with all key parameters
# - Stores optimal whitelist selection and match statistics
# - Provides configuration for subsequent pipeline steps
# 
# SAVED PARAMETERS:
# - SAMPLE: Sample identifier
# - WHITELIST: Path to optimal whitelist file
# - BEST_WHITELIST: Orientation (original/reverse_complement)
# - MATCH_RATE: Best match percentage achieved
# - BC_ERROR_THRESHOLD: Recommended error correction threshold
# - Match statistics and quality metrics
# 
# OUTPUT FORMAT:
# - Key=Value pairs for easy parsing by downstream scripts
# - Human-readable format for manual inspection
# =============================================================================
cat > "$OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt" << EOF
SAMPLE=$SAMPLE
WHITELIST=$WHITELIST
BEST_WHITELIST=$BEST_WHITELIST
MATCH_RATE=$MATCH_RATE
BC_ERROR_THRESHOLD=$BC_ERROR_THRESHOLD
MATCHES_ORIG=$MATCHES_ORIG
MATCHES_RC=$MATCHES_RC
TEST_BC_COUNT=$SMALL_BC_COUNT
SAMPLE_BARCODE=$SAMPLE_SEQ
N_COUNT=$N_COUNT
EOF

echo "Results saved to: $OUTPUT_DIR/qc/${SAMPLE}_barcode_test_results.txt"

# =============================================================================
# COMPLETION SUMMARY
# PURPOSE: Log successful completion with key results
# 
# OUTPUT: Summary of optimal whitelist selection and match rate achieved
# =============================================================================
echo "========================================="
echo "Step 2 complete for $SAMPLE"
echo "Best whitelist: $BEST_WHITELIST (${MATCH_RATE}% match rate)"
echo "End time: $(date)"
echo "=========================================="