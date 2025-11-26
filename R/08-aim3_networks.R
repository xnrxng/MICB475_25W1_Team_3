# This script computes networks under aim 3.
# Date: November 25, 2025
# Usage: Rscript R/08-aim3_networks.R

#libraries
library(here)
library(tidyverse)
library(microeco)
library(cowplot)
library(file2meco)
library(picante)
library(GUniFrac)
library(ggalluvial)
library(ggh4x)
library(ggpubr)
library(randomForest)
library(igraph)
library(rgexf)
library(htmlwidgets)
set.seed(2025)

main <- function(){
  
  # convert qiime2 data into microeco
  meco_data <- qiime2meco(
    ASV_data=here::here("moving pictures files","moving-table.qza"),
    phylo_tree=here::here("moving pictures files","moving-rooted-tree.qza"),
    taxonomy_data=here::here("moving pictures files","moving-taxonomy.qza"),
    sample_data = here::here("moving pictures files","moving-sample-metadata.tsv"))
  
}

# helper functions
calculate_relative_abundance <- function(x) x/sum(x)

main()