#!/bin/bash

# This script processed the FASTQ files using QIIME2
# Date: October 10th 2025

### exit on errors
set -euo pipefail

### activate conda env
### this script should be run in the server as it takes some absolute paths
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qiime2-amplicon-2025.4

cd data/

qiime dada2 denoise-single --i-demultiplexed-seqs data_raw/demux_seqs.qza --p-trim-left 0 --p-trunc-len 250 --o-representative-sequences data_processed/rep-seqs.qza --o-table data_processed/table.qza --o-denoising-stats data_processed/stats.qza

qiime feature-table summarize \
--i-table data_processed/table.qza \
--o-visualization data_processed/table.qzv \
--m-sample-metadata-file data_raw/colombia_metadata.txt

qiime feature-table tabulate-seqs \
--i-data data_processed/rep-seqs.qza \
--o-visualization data_processed/rep-seqs.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads data_processed/rep-seqs.qza \
  --o-classification data_processed/taxonomy.qza

qiime metadata tabulate \
  --m-input-file data_processed/taxonomy.qza \
  --o-visualization data_processed/taxonomy.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences data_processed/rep-seqs.qza \
--o-alignment data_processed/aligned-rep-seqs.qza \
--o-masked-alignment data_processed/masked-aligned-rep-seqs.qza \
--o-tree data_processed/unrooted-tree.qza \
--o-rooted-tree data_processed/rooted-tree.qza

qiime diversity alpha-rarefaction \
--i-table data_processed/table.qza \
--i-phylogeny data_processed/rooted-tree.qza \
--p-max-depth 100000 \
--m-metadata-file data_raw/colombia_metadata.txt \
--o-visualization data_processed/alpha-rarefaction.qzv

qiime tools export --input-path data_processed/table.qza --output-path data_processed/

biom convert -i data_processed/feature-table.biom --to-tsv -o data_processed/feature-table.txt

qiime tools export --input-path data_processed/taxonomy.qza --output-path data_processed/

qiime tools export --input-path data_processed/rooted-tree.qza --output-path data_processed/