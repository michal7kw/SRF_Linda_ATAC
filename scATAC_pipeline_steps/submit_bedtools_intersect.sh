#!/bin/bash

# SLURM submission script for bedtools intersect operation
# This script submits the 03_bedtools_intersect.sh job to SLURM

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps"
cd "$SCRIPT_DIR"

echo "=== Submitting Bedtools Intersect Job ==="
echo "Script directory: $SCRIPT_DIR"
echo "Timestamp: $(date)"
echo ""

# Check if the script exists
if [[ ! -f "03_bedtools_intersect.sh" ]]; then
    echo "ERROR: 03_bedtools_intersect.sh not found in $SCRIPT_DIR"
    exit 1
fi

# Make sure the script is executable
chmod +x 03_bedtools_intersect.sh

# Submit the job
echo "Submitting job with sbatch..."
job_id=$(sbatch 03_bedtools_intersect.sh | awk '{print $4}')

if [[ -n "$job_id" ]]; then
    echo "Job submitted successfully!"
    echo "Job ID: $job_id"
    echo ""
    echo "Monitor job status with:"
    echo "  squeue -j $job_id"
    echo ""
    echo "View job output with:"
    echo "  tail -f logs/03_bedtools_intersect_1.out"
    echo "  tail -f logs/03_bedtools_intersect_1.err"
    echo ""
    echo "Cancel job if needed with:"
    echo "  scancel $job_id"
else
    echo "ERROR: Failed to submit job"
    exit 1
fi

echo ""
echo "=== Job Submission Complete ==="