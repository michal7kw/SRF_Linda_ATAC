#!/bin/bash
#SBATCH --job-name=extract_barcodes
#SBATCH --output=logs/extract_barcodes_%A_%a.out
#SBATCH --error=logs/extract_barcodes_%A_%a.err
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=workq

set -euo pipefail

# Define samples
SAMPLES=("R26-Nestin-Ctrl-adult" "R26-Nestin-Mut-adult")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Extracting 16bp R2 barcode for sample: $SAMPLE"
echo "Date: $(date)"
echo "=========================================="

# Path to the barcode extraction script
EXTRACT_SCRIPT="./ATAC_data/scripts/extract_16bp_R2_barcode_fixed.sh"

# Execute the barcode extraction script
"${EXTRACT_SCRIPT}" "${SAMPLE}"

echo "=========================================="
echo "Barcode extraction completed for $SAMPLE"
echo "=========================================="