#!/bin/bash

# This script performs PICRUST.
# Date: November 19th 2025.
# Usage: bash bin/3-piicrust.sh

### exit on errors
set -euo pipefail

### activate conda env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2025.4

### filter out mitochondria and chloroplast
qiime taxa filter-table --i-table data/data_processed/table.qza --i-taxonomy data/data_processed/taxonomy.qza --p-exclude mitochondria,chloroplast --o-filtered-table data/data_processed/table-no-mitochondria-no-chloroplast.qza

### remove features with 5 or lower counts
qiime feature-table filter-features \
  --i-table data/data_processed/table-no-mitochondria-no-chloroplast.qza \
  --p-min-frequency 5 \
  --o-filtered-table data/data_processed/feature-frequency-filtered-table.qza
  
### export created files
qiime tools export \
   --input-path data/data_processed/feature-frequency-filtered-table.qza \
   --output-path data/data_processed/picrust/

qiime tools export \
   --input-path data/data_processed/rep-seqs.qza \
   --output-path data/data_processed/picrust/
   
   
### activate picrust environment
conda deactivate
conda activate picrust2

mkdir -p results/aim4/

### run picrust
picrust2_pipeline.py \
    -s data/data_processed/picrust/dna-sequences.fasta \
    -i data/data_processed/picrust/feature-table.biom \
    -o data/data_processed/picrust/team03_out
