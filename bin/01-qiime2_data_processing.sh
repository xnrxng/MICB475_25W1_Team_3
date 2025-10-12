#!/bin/bash

# This script processes the FASTQ files using QIIME2.
# Date: October 10th 2025

### exit on errors
set -euo pipefail

### activate conda env
### this script should be run in the server as it takes some absolute paths
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2025.4

### create a dir for where the data will be stored
mkdir -p data/data_raw
mkdir -p data/data_processed
cd data

### copy metadata file to data_raw
cp /datasets/project_2/colombia/colombia_metadata.txt data_raw

### import data and demultiplex
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/colombia/colombia_manifest.tsv \
  --output-path data_raw/demux_seqs.qza

### convert into qzv file
qiime demux summarize \
--i-data data_raw/demux_seqs.qza \
--o-visualization data_raw/demux.qzv
