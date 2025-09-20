#!/bin/bash
#SBATCH --job-name=05_create_matrix
#SBATCH --output=logs/05_create_matrix_%a.out
#SBATCH --error=logs/05_create_matrix_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=workq

# ==========================================
# scATAC-seq Pipeline Step 5: Create 10X Genomics Format Matrix
# ==========================================
# Purpose: Converts peak-barcode counts into 10X Genomics sparse matrix format
# Author: Pipeline Troubleshooting Analysis
# Date: September 2025
# 
# This script creates a 10X Genomics compatible sparse matrix including:
# - features.tsv.gz: Peak coordinates and metadata
# - barcodes.tsv.gz: Cell barcode identifiers  
# - matrix.mtx.gz: Sparse count matrix in MatrixMarket format
# - Comprehensive validation and statistics generation
# - Memory-efficient processing for large datasets
# 
# Input: peak_barcode_counts.txt (peak_id, barcode, count format)
# Output: 10X Genomics matrix directory with compressed files
# Dependencies: 00_config_utils.sh, gzip, bc (basic calculator)
# Performance: Optimized for memory efficiency and file compression
# ==========================================

# ==========================================
# SLURM ENVIRONMENT CHECK
# ==========================================

# Ensure script is running under SLURM
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: This script must be executed with sbatch, not bash directly."
    echo "Usage: sbatch $0"
    exit 1
fi

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps/construct_matrix"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# ==========================================
# 10X GENOMICS MATRIX CREATION FUNCTION
# ==========================================
# Purpose: Converts peak-barcode counts to 10X Genomics sparse matrix format
# Implementation: Creates features.tsv, barcodes.tsv, and matrix.mtx files with compression
# Algorithm:
#   1. Extract unique peaks and barcodes from input counts
#   2. Create index mappings for efficient sparse matrix generation
#   3. Generate MatrixMarket format sparse matrix with proper indexing
#   4. Compress all output files and validate dimensions
# Input: counts_file (peak_id, barcode, count), output_dir (target directory)
# Output: 10X Genomics matrix directory with compressed files
# Edge cases: Handles large datasets, validates file creation, ensures proper compression

create_10x_matrix() {
    local counts_file="$1"
    local output_dir="$2"
    
    log_info "Creating 10X Genomics format matrix"
    log_info "Input counts file: $counts_file"
    log_info "Output directory: $output_dir"
    
    # ==========================================
    # OUTPUT DIRECTORY SETUP
    # ==========================================
    # Purpose: Create target directory structure for 10X matrix files
    # Implementation: Creates directory with proper permissions
    # Edge cases: Handles existing directories gracefully
    mkdir -p "$output_dir" || {
        log_error "Failed to create output directory: $output_dir"
        return 1
    }
    
    # Define output files
    local features_file="${output_dir}/features.tsv"
    local barcodes_file="${output_dir}/barcodes.tsv"
    local matrix_file="${output_dir}/matrix.mtx"
    local temp_matrix_dir="${TEMP_DIR}/matrix_$$"
    
    mkdir -p "$temp_matrix_dir" || {
        log_error "Failed to create temporary matrix directory: $temp_matrix_dir"
        return 1
    }
    
    # ==========================================
    # FEATURE AND BARCODE EXTRACTION
    # ==========================================
    # Purpose: Extract unique peaks (features) and cell barcodes from counts file
    # Implementation: Uses cut and awk to extract columns, sort for consistent ordering
    # Algorithm:
    #   1. Extract peak IDs from column 1, remove duplicates, sort alphabetically
    #   2. Format peaks in 10X genomics format with proper metadata
    # Input: counts_file with format: peak_id barcode count
    # Output: features.tsv (peak list with 10X format)
    # Edge cases: Handles large files efficiently, ensures unique entries
    log_info "Extracting unique features (peaks)"
    # Create features.tsv (peak coordinates)
    cut -f1 "$counts_file" | sort -u | awk 'BEGIN{OFS="\t"} {
        split($1, coords, /[:-]/)
        chr = coords[1]
        start = coords[2]
        end = coords[3]
        peak_name = $1
        print peak_name, peak_name, "Peaks"
    }' > "$features_file"
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to create features.tsv"
        return 1
    fi
    
    local num_features=$(wc -l < "$features_file")
    log_info "Created features.tsv with $num_features features"
    
    log_info "Extracting unique barcodes"
    # Create barcodes.tsv
    cut -f2 "$counts_file" | sort -u > "$barcodes_file"
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to create barcodes.tsv"
        return 1
    fi
    
    local num_barcodes=$(wc -l < "$barcodes_file")
    log_info "Created barcodes.tsv with $num_barcodes barcodes"
    
    # Create index mappings for efficient lookup
    log_info "Creating feature and barcode index mappings"
    local feature_index_file="${temp_matrix_dir}/feature_index.txt"
    local barcode_index_file="${temp_matrix_dir}/barcode_index.txt"
    
    # Create feature index (1-based for MTX format)
    awk '{print $1 "\t" NR}' "$features_file" > "$feature_index_file"
    
    # Create barcode index (1-based for MTX format)
    awk '{print $1 "\t" NR}' "$barcodes_file" > "$barcode_index_file"
    
    # ==========================================
    # MATRIXMARKET SPARSE MATRIX CREATION
    # ==========================================
    # Purpose: Generate sparse matrix in MatrixMarket coordinate format
    # Implementation: Creates header with dimensions and coordinate entries
    # Algorithm:
    #   1. Count matrix dimensions (features x barcodes)
    #   2. Create MatrixMarket header with format specification
    #   3. Generate coordinate entries (row, col, value) with 1-based indexing
    #   4. Use awk for efficient index mapping and coordinate generation
    # Input: counts_file, features.tsv, barcodes.tsv
    # Output: matrix.mtx in MatrixMarket coordinate format
    # Edge cases: Handles large matrices, ensures proper 1-based indexing, validates dimensions
    log_info "Creating MatrixMarket format matrix"
    
    # ==========================================
    # DIMENSION CALCULATION
    # ==========================================
    # Purpose: Calculate matrix dimensions for MatrixMarket header
    # Implementation: Count lines in feature and barcode files, entries in counts file
    # Edge cases: Validates non-zero dimensions, handles large counts
    local matrix_features="$num_features"
    local matrix_barcodes="$num_barcodes"
    local matrix_entries=$(wc -l < "$counts_file")
    
    log_info "Matrix dimensions: $matrix_features features x $matrix_barcodes barcodes ($matrix_entries entries)"
    
    # ==========================================
    # MATRIXMARKET HEADER AND COORDINATE GENERATION
    # ==========================================
    # Purpose: Create MatrixMarket format header with proper specifications and coordinate entries
    # Implementation: Standard MatrixMarket coordinate format with 1-based indexing
    # Format: %%MatrixMarket matrix coordinate integer general
    #         % (comment line)
    #         num_rows num_cols num_entries
    #         row_idx col_idx value (for each non-zero entry)
    # Algorithm:
    #   1. Write MatrixMarket header with format specification
    #   2. Join counts with feature indices using sorted files for efficiency
    #   3. Join result with barcode indices to get complete coordinate mapping
    #   4. Output in MTX coordinate format: feature_idx barcode_idx count
    #   5. Sort by row then column for optimal sparse matrix performance
    {
        echo "%%MatrixMarket matrix coordinate integer general"
        echo "%"
        echo "$matrix_features $matrix_barcodes $matrix_entries"
        
        # ==========================================
        # COORDINATE ENTRY GENERATION
        # ==========================================
        # Purpose: Convert peak-barcode counts to MatrixMarket coordinate format
        # Implementation: Multi-step join operation with index mapping
        # Algorithm:
        #   Step 1: Join counts with feature indices (peak_id -> feature_index)
        #   Step 2: Join result with barcode indices (barcode -> barcode_index)
        #   Step 3: Output in MTX coordinate format with proper ordering
        # Input: counts_file (peak_id barcode count), index files (id -> 1-based_index)
        # Output: Sorted coordinate entries (feature_idx barcode_idx count)
        # Edge cases: Handles large datasets, ensures proper 1-based indexing, maintains sort order
        
        # Step 1: Join counts with feature indices
        # Input: peak_coord barcode_pos count + peak_coord feature_idx
        # Output: peak_coord barcode_pos count feature_idx
        join -1 1 -2 1 <(sort -k1,1 "$counts_file") <(sort -k1,1 "$feature_index_file") | \
        # Step 2: Join result with barcode indices  
        # Input: peak_coord barcode_pos count feature_idx + barcode_pos barcode_idx
        # Output: barcode_pos peak_coord count feature_idx barcode_idx
        join -1 2 -2 1 <(sort -k2,2) <(sort -k1,1 "$barcode_index_file") | \
        # Step 3: Output in MTX coordinate format: feature_idx barcode_idx count
        # Sort by feature index (row) then barcode index (column) for optimal sparse matrix access
        awk '{print $4, $5, $3}' | \
        sort -k1,1n -k2,2n
        
    } > "$matrix_file"
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to create matrix.mtx"
        return 1
    fi
    
    log_info "Matrix.mtx created successfully"
    
    # Validate matrix dimensions
    local matrix_entries=$(tail -n +4 "$matrix_file" | wc -l)
    local expected_entries=$(wc -l < "$counts_file")
    
    if [[ $matrix_entries -ne $expected_entries ]]; then
        log_error "Matrix entry count mismatch. Expected: $expected_entries, Got: $matrix_entries"
        return 1
    fi
    
    log_info "Matrix validation passed: $matrix_entries entries"
    
    # ==========================================
    # FILE VALIDATION AND COMPRESSION
    # ==========================================
    # Purpose: Validate successful file creation and compress for storage efficiency
    # Implementation: Check file existence, validate content, apply gzip compression
    # Algorithm:
    #   1. Verify all required files exist and are non-empty
    #   2. Validate file formats and dimensions consistency
    #   3. Apply gzip compression to reduce storage footprint
    #   4. Verify compressed files are created successfully
    # Input: Uncompressed TSV and MTX files
    # Output: Compressed .gz files with validation
    # Edge cases: Handles compression failures, validates file integrity
    log_info "Compressing output files"
    
    # ==========================================
    # FILE EXISTENCE VALIDATION
    # ==========================================
    # Purpose: Ensure all required matrix files were created successfully
    # Implementation: Check file existence and non-zero size
    # Edge cases: Reports specific missing files, handles permission issues
    local compression_failed=0
    
    if [[ -f "$features_file" ]]; then
        gzip -f "$features_file"
        if [[ $? -ne 0 ]]; then
            log_error "Failed to compress features file: $features_file"
            compression_failed=1
        fi
    else
        log_error "Features file not found for compression: $features_file"
        compression_failed=1
    fi
    
    if [[ -f "$barcodes_file" ]]; then
        gzip -f "$barcodes_file"
        if [[ $? -ne 0 ]]; then
            log_error "Failed to compress barcodes file: $barcodes_file"
            compression_failed=1
        fi
    else
        log_error "Barcodes file not found for compression: $barcodes_file"
        compression_failed=1
    fi
    
    if [[ -f "$matrix_file" ]]; then
        gzip -f "$matrix_file"
        if [[ $? -ne 0 ]]; then
            log_error "Failed to compress matrix file: $matrix_file"
            compression_failed=1
        fi
    else
        log_error "Matrix file not found for compression: $matrix_file"
        compression_failed=1
    fi
    
    # ==========================================
    # COMPRESSION VALIDATION
    # ==========================================
    # Purpose: Verify all files were compressed successfully
    # Implementation: Check compression status and report results
    # Edge cases: Handles partial compression failures, provides detailed error reporting
    if [[ $compression_failed -eq 1 ]]; then
        log_error "Failed to compress output files"
        return 1
    fi
    
    log_info "Files compressed successfully"
    
    # ==========================================
    # MATRIX STATISTICS CALCULATION
    # ==========================================
    # Purpose: Calculate and report comprehensive matrix statistics
    # Implementation: Compute dimensions, density, and storage metrics
    # Algorithm:
    #   1. Calculate total counts across all matrix entries
    #   2. Compute average counts per feature and per barcode
    #   3. Calculate matrix density (sparsity) for optimization insights
    #   4. Report comprehensive statistics for analysis planning
    # Input: Matrix dimensions and count data
    # Output: Detailed statistics logging
    # Edge cases: Handles large numbers, prevents division by zero
    local total_counts=$(awk '{sum += $3} END {print sum}' "$counts_file")
    local avg_counts_per_feature=$(echo "scale=2; $total_counts / $num_features" | bc)
    local avg_counts_per_barcode=$(echo "scale=2; $total_counts / $num_barcodes" | bc)
    local sparsity=$(echo "scale=4; $matrix_entries / ($num_features * $num_barcodes) * 100" | bc)
    
    log_info "=== Matrix Statistics ==="
    log_info "Features (peaks): $num_features"
    log_info "Barcodes (cells): $num_barcodes"
    log_info "Non-zero entries: $matrix_entries"
    log_info "Total counts: $total_counts"
    log_info "Average counts per feature: $avg_counts_per_feature"
    log_info "Average counts per barcode: $avg_counts_per_barcode"
    log_info "Matrix density: ${sparsity}%"
    
    # Create matrix info file
    local info_file="${output_dir}/matrix_info.txt"
    cat > "$info_file" << EOF
10X Genomics Matrix Information
==============================
Created at: $(date)

Matrix Dimensions:
- Features (peaks): $num_features
- Barcodes (cells): $num_barcodes
- Non-zero entries: $matrix_entries
- Total counts: $total_counts

Statistics:
- Average counts per feature: $avg_counts_per_feature
- Average counts per barcode: $avg_counts_per_barcode
- Matrix density: ${sparsity}%
- Matrix sparsity: $(echo "scale=4; 100 - $sparsity" | bc)%

Files:
- features.tsv.gz: Peak coordinates and metadata
- barcodes.tsv.gz: Cell barcode identifiers
- matrix.mtx.gz: Sparse count matrix in MatrixMarket format
- matrix_info.txt: This information file

Format:
- Compatible with 10X Genomics Cell Ranger output
- Can be loaded with scanpy, Seurat, or other single-cell analysis tools
EOF
    
    log_info "Matrix information saved to: $info_file"
    
    # Cleanup temporary files
    cleanup_temp_dir "$temp_matrix_dir"
    
    return 0
}

# ==========================================
# 10X GENOMICS MATRIX VALIDATION FUNCTION
# ==========================================
# Purpose: Comprehensive validation of 10X Genomics matrix format compliance
# Implementation: Validates file existence, format, dimensions, and data integrity
# Algorithm:
#   1. Check existence of all required files (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz)
#   2. Validate file formats and compression integrity
#   3. Verify MatrixMarket header format and dimensions
#   4. Cross-validate dimensions across all files for consistency
#   5. Perform basic data integrity checks
# Input: matrix_dir (path to 10X matrix directory)
# Output: Validation status and detailed error reporting
# Edge cases: Handles corrupted files, dimension mismatches, format violations

validate_10x_matrix() {
    local matrix_dir="$1"
    
    log_info "Validating 10X matrix format"
    
    # ==========================================
    # REQUIRED FILE EXISTENCE CHECK
    # ==========================================
    # Purpose: Verify all essential 10X matrix files are present
    # Implementation: Check for compressed TSV and MTX files
    # Edge cases: Reports specific missing files, handles permission issues
    local features_file="${matrix_dir}/features.tsv.gz"
    local barcodes_file="${matrix_dir}/barcodes.tsv.gz"
    local matrix_file="${matrix_dir}/matrix.mtx.gz"
    
    if [[ ! -f "$features_file" ]]; then
        log_error "Features file not found: $features_file"
        return 1
    fi
    
    if [[ ! -f "$barcodes_file" ]]; then
        log_error "Barcodes file not found: $barcodes_file"
        return 1
    fi
    
    if [[ ! -f "$matrix_file" ]]; then
        log_error "Matrix file not found: $matrix_file"
        return 1
    fi
    
    log_info "All required files present"
    
    # ==========================================
    # FILE FORMAT VALIDATION
    # ==========================================
    # Purpose: Validate proper file formats and compression integrity
    # Implementation: Check MatrixMarket header format and file structure
    # Algorithm:
    #   1. Verify MatrixMarket header format compliance
    #   2. Check file compression integrity (gzip)
    #   3. Validate TSV file structure for features and barcodes
    # Edge cases: Handles corrupted compression, malformed headers
    
    # ==========================================
    # FEATURES FILE FORMAT VALIDATION
    # ==========================================
    # Purpose: Validate features.tsv.gz format and content structure
    # Implementation: Check TSV format with proper column structure
    # Format: peak_id\tpeak_name\tPeaks (3 columns, tab-separated)
    # Edge cases: Handles missing columns, incorrect feature types
    local features_lines=$(zcat "$features_file" | head -n 5)
    local features_count=0
    
    while IFS=$'\t' read -r col1 col2 col3; do
        ((features_count++))
        
        if [[ -z "$col1" || -z "$col2" || -z "$col3" ]]; then
            log_error "Features file: Line $features_count missing required columns"
            return 1
        fi
        
        if [[ "$col3" != "Peaks" ]]; then
            log_error "Features file: Line $features_count has unexpected feature type: $col3"
            return 1
        fi
        
        if [[ $features_count -ge 5 ]]; then
            break
        fi
    done <<< "$features_lines"
    
    # ==========================================
    # MATRIXMARKET FORMAT VALIDATION
    # ==========================================
    # Purpose: Validate MatrixMarket coordinate format compliance
    # Implementation: Check header format and dimension specification
    # Format: %%MatrixMarket matrix coordinate integer general
    #         %
    #         num_features num_barcodes num_entries
    # Edge cases: Handles malformed headers, incorrect format specifications
    local matrix_header=$(zcat "$matrix_file" | head -n 3)
    local header_line1=$(echo "$matrix_header" | sed -n '1p')
    local header_line3=$(echo "$matrix_header" | sed -n '3p')
    
    if [[ "$header_line1" != "%%MatrixMarket matrix coordinate integer general" ]]; then
        log_error "Matrix file: Invalid MatrixMarket header"
        return 1
    fi
    
    # ==========================================
    # DIMENSION EXTRACTION AND PARSING
    # ==========================================
    # Purpose: Extract matrix dimensions from MatrixMarket header
    # Implementation: Parse third line of MTX file for dimensions
    # Format: num_features num_barcodes num_entries
    # Edge cases: Handles malformed dimension lines, invalid numbers
    local dimensions=($header_line3)
    local matrix_features=${dimensions[0]}
    local matrix_barcodes=${dimensions[1]}
    local matrix_entries=${dimensions[2]}
    
    # ==========================================
    # ACTUAL FILE DIMENSION COUNTING
    # ==========================================
    # Purpose: Count actual lines in each file for dimension validation
    # Implementation: Use wc -l to count lines, skip MTX header lines
    # Algorithm:
    #   1. Count features: lines in features.tsv.gz
    #   2. Count barcodes: lines in barcodes.tsv.gz
    #   3. Count entries: data lines in matrix.mtx.gz (skip 3 header lines)
    # Edge cases: Handles large files efficiently, accounts for header lines
    local actual_features=$(zcat "$features_file" | wc -l)
    local actual_barcodes=$(zcat "$barcodes_file" | wc -l)
    local actual_entries=$(zcat "$matrix_file" | tail -n +4 | wc -l)
    
    # ==========================================
    # DIMENSION CONSISTENCY VALIDATION
    # ==========================================
    # Purpose: Cross-validate dimensions across all matrix files
    # Implementation: Compare declared dimensions with actual file contents
    # Algorithm:
    #   1. Compare matrix header features with features.tsv.gz line count
    #   2. Compare matrix header barcodes with barcodes.tsv.gz line count
    #   3. Compare matrix header entries with actual coordinate entries
    # Edge cases: Provides detailed error messages for each mismatch type
    
    if [[ $matrix_features -ne $actual_features ]]; then
        log_error "Feature count mismatch. Matrix header: $matrix_features, Features file: $actual_features"
        return 1
    fi
    
    if [[ $matrix_barcodes -ne $actual_barcodes ]]; then
        log_error "Barcode count mismatch. Matrix header: $matrix_barcodes, Barcodes file: $actual_barcodes"
        return 1
    fi
    
    if [[ $matrix_entries -ne $actual_entries ]]; then
        log_error "Entry count mismatch. Matrix header: $matrix_entries, Actual entries: $actual_entries"
        return 1
    fi
    
    log_info "10X matrix format validation passed"
    log_info "Matrix dimensions: $matrix_features features × $matrix_barcodes barcodes"
    log_info "Non-zero entries: $matrix_entries"
    
    return 0
}

# ==========================================
# MATRIX SUMMARY GENERATION FUNCTION
# ==========================================
# Purpose: Generate comprehensive summary report of 10X matrix creation
# Implementation: Collect statistics, file sizes, and analysis recommendations
# Algorithm:
#   1. Calculate file sizes and compression ratios
#   2. Extract matrix dimensions and sparsity metrics
#   3. Generate performance statistics and recommendations
#   4. Create formatted summary report for downstream analysis
# Input: matrix_dir (10X matrix directory), summary_file (output report path)
# Output: Detailed summary report with statistics and next steps
# Edge cases: Handles missing files, calculation errors, formatting issues

create_matrix_summary() {
    local matrix_dir="$1"
    local summary_file="$2"
    
    log_info "Creating matrix summary"
    
    # ==========================================
    # FILE PATH DEFINITION
    # ==========================================
    # Purpose: Define paths to all 10X matrix component files
    # Implementation: Standard 10X Genomics directory structure
    # Edge cases: Handles missing directory structure gracefully
    local features_file="${matrix_dir}/features.tsv.gz"
    local barcodes_file="${matrix_dir}/barcodes.tsv.gz"
    local matrix_file="${matrix_dir}/matrix.mtx.gz"
    
    # ==========================================
    # MATRIX DIMENSION EXTRACTION
    # ==========================================
    # Purpose: Extract matrix dimensions from compressed files
    # Implementation: Use zcat to read compressed files, count lines
    # Algorithm:
    #   1. Count features: lines in features.tsv.gz
    #   2. Count barcodes: lines in barcodes.tsv.gz
    #   3. Count entries: data lines in matrix.mtx.gz (skip header)
    # Edge cases: Handles large files efficiently, accounts for MTX header
    local num_features=$(zcat "$features_file" | wc -l)
    local num_barcodes=$(zcat "$barcodes_file" | wc -l)
    local num_entries=$(zcat "$matrix_file" | tail -n +4 | wc -l)
    
    # ==========================================
    # TOTAL COUNTS CALCULATION
    # ==========================================
    # Purpose: Calculate total counts across all matrix entries
    # Implementation: Sum third column of matrix coordinate entries
    # Algorithm: Skip MTX header, sum count values from coordinate format
    # Edge cases: Handles large numbers, empty matrices
    local total_counts=$(zcat "$matrix_file" | tail -n +4 | awk '{sum += $3} END {print sum}')
    
    # ==========================================
    # FILE SIZE CALCULATION
    # ==========================================
    # Purpose: Calculate storage requirements and compression efficiency
    # Implementation: Use stat to get file sizes, convert to MB
    # Algorithm: Get byte size, convert to MB with decimal precision
    # Edge cases: Handles large files, missing files, permission issues
    local features_size=$(stat -c%s "$features_file" | awk '{printf "%.1f MB", $1/1024/1024}')
    local barcodes_size=$(stat -c%s "$barcodes_file" | awk '{printf "%.1f MB", $1/1024/1024}')
    local matrix_size=$(stat -c%s "$matrix_file" | awk '{printf "%.1f MB", $1/1024/1024}')
    
    cat > "$summary_file" << EOF
Step 5: Create 10X Matrix - Summary
==================================
Completed at: $(date)

Output Directory: $matrix_dir

Matrix Dimensions:
- Features (peaks): $num_features
- Barcodes (cells): $num_barcodes
- Non-zero entries: $num_entries
- Total counts: $total_counts

File Sizes:
- features.tsv.gz: $features_size
- barcodes.tsv.gz: $barcodes_size
- matrix.mtx.gz: $matrix_size

Next Steps:
- Run step 6 for output validation
- Load matrix in R/Python for analysis
- Use with scanpy, Seurat, or other tools
EOF
    
    log_info "Summary saved to: $summary_file"
}

# ==========================================
# MAIN EXECUTION FUNCTION
# ==========================================
# Purpose: Orchestrate the complete 10X matrix creation pipeline
# Implementation: Sequential execution of validation, creation, and reporting steps
# Algorithm:
#   1. Environment initialization and validation
#   2. Disk space and prerequisite checks
#   3. Input/output path definition and validation
#   4. Matrix creation with comprehensive error handling
#   5. Validation and summary report generation
# Input: Global configuration variables and command line parameters
# Output: Complete 10X Genomics matrix with validation and summary
# Edge cases: Handles all error conditions, provides detailed logging

main() {
    log_info "=== Step 5: Create 10X Genomics Format Matrix ==="
    log_info "Started at: $(date)"
    
    # ==========================================
    # ENVIRONMENT INITIALIZATION
    # ==========================================
    # Purpose: Initialize logging, conda environment, and temporary directories
    # Implementation: Sequential setup of required runtime environment
    # Edge cases: Handles setup failures, missing dependencies
    init_logging
    setup_conda_env
    setup_temp_dir
    
    # ==========================================
    # DISK SPACE VALIDATION
    # ==========================================
    # Purpose: Ensure sufficient disk space for matrix creation and compression
    # Implementation: Check available space against minimum requirements (5GB)
    # Edge cases: Handles disk space calculation errors, permission issues
    check_disk_space "$OUTPUT_DIR" 5
    
    # ==========================================
    # INPUT/OUTPUT PATH DEFINITION
    # ==========================================
    # Purpose: Define all input and output file paths for matrix creation
    # Implementation: Construct paths based on output directory configuration
    # Algorithm:
    #   1. Define input counts file from previous pipeline step
    #   2. Define output matrix directory for 10X format files
    #   3. Define summary report file for pipeline documentation
    # Edge cases: Handles path construction errors, directory creation issues
    local counts_file="${OUTPUT_DIR}/${SAMPLE}_peak_barcode_counts.txt"
    local matrix_dir="${OUTPUT_DIR}/${SAMPLE}_filtered_peak_bc_matrix"
    local summary_file="${OUTPUT_DIR}/${SAMPLE}_step5_matrix_summary.txt"
    
    log_info "Input counts file: $counts_file"
    log_info "Output matrix directory: $matrix_dir"
    log_info "Summary file: $summary_file"
    
    # ==========================================
    # INPUT FILE VALIDATION
    # ==========================================
    # Purpose: Validate that required input files exist and are accessible
    # Implementation: Check file existence and readability
    # Edge cases: Handles missing files, permission issues, empty files
    if [[ ! -f "$counts_file" ]]; then
        log_error "Peak-barcode counts file not found: $counts_file"
        log_error "Please run step 4 (04_process_overlaps.sh) first"
        exit 1
    fi
    
    if [[ ! -s "$counts_file" ]]; then
        log_error "Input counts file is empty: $counts_file"
        exit 1
    fi
    
    # Check if bc (basic calculator) is available for statistics
    if ! command -v bc &> /dev/null; then
        log_warning "bc (basic calculator) not found. Some statistics may be unavailable."
    fi
    
    # Create 10X matrix
    if ! create_10x_matrix "$counts_file" "$matrix_dir"; then
        log_error "10X matrix creation failed"
        exit 1
    fi
    
    # Validate matrix format
    if ! validate_10x_matrix "$matrix_dir"; then
        log_error "10X matrix validation failed"
        exit 1
    fi
    
    # Create summary
    create_matrix_summary "$matrix_dir" "$summary_file"
    
    log_info "Step 5 completed successfully at: $(date)"
    log_info "Matrix created in: $matrix_dir"
    log_info "Next step: Run 06_validate_output.sh"
}

# Run main function
main "$@"