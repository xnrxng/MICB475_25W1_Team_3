#!/bin/bash

# This script downloads the FASTQ files from SRA.
# Date: November 3rd 2025
# Usage: bash bin/0-data_download.sh

### exit on errors
set -euo pipefail

### activate conda env
### this script should be run in the server as it takes some absolute paths
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2025.4

### for loop that downloads each fastq file
input="data/data_raw/srr_ids.txt"
outdir="/datasets/project_2/colombia/seqs/"

mkdir -p "$outdir"

while read srr; do
  echo "Processing $srr ..."

  prefix=${srr:0:6}                    
  subdir=$(printf "%03d" ${srr: -3})   
  base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${subdir}/${srr}"

  wget -c "${base}/${srr}_1.fastq.gz" -P "$outdir" || echo "No _1 file for $srr"

done < "$input"

echo "download complete. "