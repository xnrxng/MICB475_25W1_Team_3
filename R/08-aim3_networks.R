# This script computes networks under aim 3.
# Date: November 25, 2025
# Usage: Rscript R/08-aim3_networks.R

#libraries
library(tidyverse)
library(microeco)
library(cowplot)
library(file2meco)
library(picante)
library(GUniFrac)
library(ggalluvial)
library(ggh4x)
library(ggpubr)
library(igraph)
library(htmlwidgets)
library(dplyr)
library(rgexf)
set.seed(2025)

main <- function(){
  
  # convert qiime2 data into microeco
  meco_data <- qiime2meco(
    feature_table="data/data_processed/table-no-mitochondria-no-chloroplast.qza",
    phylo_tree="data/data_processed/rooted-tree.qza",
    taxonomy_table ="data/data_processed/taxonomy.qza",
    sample_table = "data/data_raw/colombia_metadata.txt")
  
  # update metadata to include groups
  meta <- meco_data$sample_table
  
  # update / append new grouping variables
  meta <- meta |>
    mutate(
      fibre_group = ifelse(fiber >= 20, "high", "low"),
      exercise_group = ifelse(MET_mins_per_week > 1000, "high", "low"),
      fibre_exercise_group = paste(fibre_group, exercise_group, sep = "_")
    )
  
  # write it back to the microeco object
  meco_data$sample_table <- meta
  
  #filter
  #Create a dummy copy of our dataset
  dataset_healthy <- clone(meco_data)
  #Filter the samples to only include healthy samples
  dataset_healthy$sample_table <- dataset_healthy$sample_table |> filter( Cardiometabolic_status== "Healthy")
  #tidy the whole microtable
  dataset_healthy$tidy_dataset()
  
  #Create a dummy copy of our dataset
  dataset_abnormal <- clone(meco_data)
  #Filter the samples to only include abnormal samples
  dataset_abnormal$sample_table <- dataset_abnormal$sample_table |> filter( Cardiometabolic_status== "Abnormal")
  #tidy the whole microtable
  dataset_abnormal$tidy_dataset()
  
  # network
  t_healthy <- trans_network$new(dataset = dataset_healthy,
                          cal_cor = "base",
                          taxa_level = "Genus",
                          filter_thres = 0.0001,
                          cor_method = "spearman")
  
  t_abnormal <- trans_network$new(dataset = dataset_abnormal,
                                 cal_cor = "base",
                                 taxa_level = "Genus",
                                 filter_thres = 0.0001,
                                 cor_method = "spearman")
  
  # construct the actual network and add modules
  t_healthy$cal_network(p_thres = 0.05, COR_optimization = TRUE)
  t_healthy$cal_module()
  
  t_abnormal$cal_network(p_thres = 0.05, COR_optimization = TRUE)
  t_abnormal$cal_module()
  
  t_healthy$save_network("results/aim3/networks/01-healthy_network.gexf")
  t_abnormal$save_network("results/aim3/networks/02-abnormal_network.gexf")
  
  #calculate network attributes
  t_healthy$cal_network_attr()
  t_abnormal$cal_network_attr()
  
  # classify the node
  t_healthy$get_node_table(node_roles = TRUE)
  t_healthy$get_edge_table()
  t_healthy$get_adjacency_matrix()
  
  t_abnormal$get_node_table(node_roles = TRUE)
  t_abnormal$get_edge_table()
  t_abnormal$get_adjacency_matrix()
  
  #plot roles 
  role_cols <- c(
    "Peripheral nodes" = "#1f77b4",
    "Connectors"  = "#ff7f0e",
    "Module hubs" = "#2ca02c",
    "NA" = "lightgray"
  )
  
  t_healthy$plot_taxa_roles(use_type = 1)+
    scale_color_manual(values = role_cols)
  ggsave("results/aim3/networks/07-healthy_node_roles.png")
  
  t_abnormal$plot_taxa_roles(use_type = 1)+
    scale_color_manual(values = role_cols)
  ggsave("results/aim3/networks/08-abnormal_node_roles.png")
  
  #plot roles based on phylum
  role_shapes <- c(
    "Peripheral nodes" = 17,   
    "Connectors"  = 16,   
    "Module hubs" = 15,   
    "NA" = 18   
  )
  
  t_healthy$plot_taxa_roles(use_type = 2) +
    scale_shape_manual(values = role_shapes, name = "Taxa role") +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.63, size = 10))
  
  ggsave("results/aim3/networks/09-healthy_node_roles_phylum.png", width = 8, height = 5)
  
  t_abnormal$plot_taxa_roles(use_type = 2) +
    scale_shape_manual(values = role_shapes, name = "Taxa role") +
    theme_bw() +
    ylab(NULL) +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.63, size = 10))
  
  ggsave("results/aim3/networks/10-abnormal_node_roles_phylum.png", width = 8, height = 5)
  
  #eigengene analysis of modules
  t_healthy$cal_eigen()
  t_abnormal$cal_eigen()
  
  #correlation
  t2_healthy <- trans_env$new(dataset = dataset_healthy, 
                      env_cols = c("fiber", "MET_mins_per_week", "age_years", "Calorie_intake", "BMI", "adiponectin"))
  t2_healthy$cal_cor(add_abund_table = t_healthy$res_eigen)
  
  t2_healthy$plot_cor(text_y_order = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9") )
  ggsave("results/aim3/networks/11-healthy_meta_correlation.png", width = 6, height = 5)
  
  t2_abnormal <- trans_env$new(dataset = dataset_abnormal, 
                              env_cols = c("fiber", "MET_mins_per_week", "age_years", "Calorie_intake", "BMI", "adiponectin"))
  t2_abnormal$cal_cor(add_abund_table = t_abnormal$res_eigen)
  
  t2_abnormal$plot_cor(text_y_order = c("M1", "M2", "M3", "M4", "M5", "M6", "M7"),
                       text_x_order = c("BMI", "MET_mins_per_week", "adiponectin", "Calorie_intake", "fiber"), sig_label_size = 8)
  ggsave("results/aim3/networks/12-abnormal_meta_correlation.png", width = 6, height = 5)
  
  write_tsv(t_healthy[["res_node_table"]], "results/aim3/networks/13-healthy_network_roles.tsv")
  write_tsv(t_abnormal[["res_node_table"]], "results/aim3/networks/14-abnormal_network_roles.tsv")
}

# helper functions
calculate_relative_abundance <- function(x) x/sum(x)

main()