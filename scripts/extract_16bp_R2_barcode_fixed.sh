#!/bin/bash

# Check if sample name is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <sample_name>"
    exit 1
fi

SAMPLE_NAME=$1
SOURCE_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin"
OUTPUT_DATA_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/ATAC_data/DATA/nestin_processed"
mkdir -p "${OUTPUT_DATA_DIR}"

INPUT_R2_FILE="${SOURCE_DATA_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"
OUTPUT_R2_FILE="${OUTPUT_DATA_DIR}/${SAMPLE_NAME}_R2_16bp_fixed.fastq.gz"

echo "Extracting last 16bp barcode from ${INPUT_R2_FILE}..."

pigz -dc "${INPUT_R2_FILE}" | \
    awk '{
        if(NR%4==1) {print $0;}  # header
        else if(NR%4==2) {print substr($0, length($0)-15);}  # extract last 16 bases
        else if(NR%4==3) {print $0;}  # plus line
        else if(NR%4==0) {print substr($0, length($0)-15);}  # extract last 16 quality scores
    }' | pigz -c > "${OUTPUT_R2_FILE}"

echo "Generated ${OUTPUT_R2_FILE}"