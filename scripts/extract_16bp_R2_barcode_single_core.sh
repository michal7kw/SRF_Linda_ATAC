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
OUTPUT_R2_FILE="${OUTPUT_DATA_DIR}/${SAMPLE_NAME}_R2_16bp.fastq.gz"

echo "Extracting 16bp barcode from ${INPUT_R2_FILE}..."

zcat "${INPUT_R2_FILE}" | \
    awk '{
        if(NR%4==1) print $0;  # header
        else if(NR%4==2) print substr($0,7,16);  # extract bases 7-22
        else if(NR%4==3) print $0;  # plus line
        else if(NR%4==0) print substr($0,7,16);  # extract quality 7-22
    }' | gzip > "${OUTPUT_R2_FILE}"

echo "Generated ${OUTPUT_R2_FILE}"