#!/bin/bash

# This script processed the FASTQ files using QIIME2
# Date: October 10th 2025

### activate conda env
### this script should be run in the server as it takes some absolute paths
conda activate qiime2-2025.4

### create a dir for where the data will be stored
mkdir -r data/data_raw
mkdir -r data/data_processed
cd data

# import data and demultiplex
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/colombia/colombia_manifest.tsv \
  --output-path data_raw/demux_seqs.qza

qiime demux summarize \
--i-data demux_seqs.qza \
--o-visualization data_raw/demux.qzv