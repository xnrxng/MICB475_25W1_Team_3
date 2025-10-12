#!/bin/bash

# This script processed the FASTQ files using QIIME2
# Date: October 10th 2025

### activate conda env
### this script should be run in the server as it takes some absolute paths
conda activate qiime2-2025.4

qiime dada2 denoise-single --i-demultiplexed-seqs demux_seqs.qza --p-trim-left 0 --p-trunc-len 147 --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 8000 \
--m-metadata-file /mnt/datasets/project_1/moving_pictures/sample-
metadata.tsv \
--o-visualization alpha-rarefaction.qzv

qiime tools export
--input-path FILE.qza
--output-path FILE_EXPORT_DIR

biom convert
-i feature-table.biom
--to-tsv -o feature-table.txt