#!/bin/bash

# Test script to verify memory detection functionality
# This script tests the SLURM memory detection without running the full pipeline

SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/scATAC_pipeline_steps"
cd "$SCRIPT_DIR"
source "00_config_utils.sh"

# Source the memory detection function from the main script
source <(grep -A 25 "get_slurm_memory_gb()" 03_bedtools_intersect.sh)

echo "=== Memory Detection Test ==="
echo "Current SLURM environment variables:"
echo "SLURM_MEM_PER_NODE: ${SLURM_MEM_PER_NODE:-'not set'}"
echo "SLURM_MEM_PER_CPU: ${SLURM_MEM_PER_CPU:-'not set'}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-'not set'}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-'not set'}"
echo ""

# Test memory detection
echo "Testing memory detection..."
detected_mem=$(get_slurm_memory_gb)
echo "Detected SLURM memory: ${detected_mem}GB"

# Test memory calculation logic
available_mem_gb=$(free -g | awk '/^Mem:/{print int($7*0.8)}')
usable_mem_gb=$((detected_mem < available_mem_gb ? detected_mem : available_mem_gb))
sort_mem_gb=$((usable_mem_gb * 60 / 100))

echo ""
echo "Memory calculation results:"
echo "Available system memory (80%): ${available_mem_gb}GB"
echo "Usable memory for job: ${usable_mem_gb}GB"
echo "Sort memory allocation (60%): ${sort_mem_gb}GB"

if [[ $sort_mem_gb -lt 2 ]]; then
    echo "WARNING: Very low memory detected. Would use minimum 2GB."
elif [[ $sort_mem_gb -lt 4 ]]; then
    echo "WARNING: Low memory detected. Would use conservative 2GB."
else
    echo "Memory allocation looks good for processing."
fi

echo ""
echo "=== Test Complete ==="