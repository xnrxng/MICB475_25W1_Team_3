# This script performs beta diversity under aim 2. 
# Date: November 25, 2025
# Usage: Rscript R/09-aim4_picrust2.R

# Load in libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpicrust2)
library(pheatmap)
# Set the seed for reproducibility 
set.seed(2025)

main <- function(){
  #load data and get meta. prepare comparison groups
  phylo_obj <- readRDS("data/data_raw/phyloseq_object.rds")
  
  meta <- sample_data(phylo_obj) |> data.frame() |> rownames_to_column("sample_name") |>
    mutate(fibre_group = ifelse(fiber >= 20, "high", "low"),
           exercise_group = ifelse(MET_mins_per_week > 1000, "high", "low"))
  
  meta$fibre_group <- as.factor(meta$fibre_group)
  meta$exercise_group <- as.factor(meta$exercise_group)
  
  meta$fibre_group <- relevel(meta$fibre_group, ref = "low")
  meta$exercise_group <- relevel(meta$exercise_group, ref = "low")
  
  meta$Cardiometabolic_status <- as.factor(meta$Cardiometabolic_status)
  meta$Cardiometabolic_status <- relevel(meta$Cardiometabolic_status, ref = "Healthy")
  
  # read ko
  ko <- read.delim("data/data_processed/picrust/team03_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", row.names = 1)
  
  ####################################### DAA
  #healthy
  meta_healthy <- meta[meta$Cardiometabolic_status == "Healthy", ]
  ko_healthy <- ko[, colnames(ko) %in% meta_healthy$sample_name]
  
  daa_healthy_fibre <- pathway_daa(abundance = ko_healthy,
                                   metadata = meta_healthy,
                                   group = "fibre_group",
                                   daa_method = "edgeR")
  
  daa_healthy_fibre_annot <- pathway_annotation(pathway = "KO", 
                                                daa_results_df = daa_healthy_fibre,
                                                ko_to_kegg = TRUE)
  
  daa_healthy_exercise <- pathway_daa(abundance = ko_healthy,
                                   metadata = meta_healthy,
                                   group = "exercise_group",
                                   daa_method = "edgeR")
  
  daa_healthy_exercise_annot <- pathway_annotation(pathway = "KO", 
                                                daa_results_df = daa_healthy_exercise,
                                                ko_to_kegg = TRUE)

  saveRDS(daa_healthy_fibre_annot, "results/aim4/01-daa_healthy_fibre.rds")
  saveRDS(daa_healthy_exercise_annot, "results/aim4/02-daa_healthy_exercise.rds")
  
  # abnormal
  meta_abnormal <- meta[meta$Cardiometabolic_status == "Abnormal", ]
  ko_abnormal <- ko[, colnames(ko) %in% meta_abnormal$sample_name]
  
  daa_abnormal_fibre <- pathway_daa(abundance = ko_abnormal,
                                   metadata = meta_abnormal,
                                   group = "fibre_group",
                                   daa_method = "edgeR")
  
  daa_abnormal_fibre_annot <- pathway_annotation(pathway = "KO", 
                                                daa_results_df = daa_abnormal_fibre,
                                                ko_to_kegg = TRUE)
  
  daa_abnormal_exercise <- pathway_daa(abundance = ko_abnormal,
                                      metadata = meta_abnormal,
                                      group = "exercise_group",
                                      daa_method = "edgeR")
  
  daa_abnormal_exercise_annot <- pathway_annotation(pathway = "KO", 
                                                   daa_results_df = daa_abnormal_exercise,
                                                   ko_to_kegg = TRUE)

  saveRDS(daa_abnormal_fibre_annot, "results/aim4/03-daa_abnormal_fibre.rds")
  saveRDS(daa_abnormal_exercise_annot, "results/aim4/04-daa_abnormal_exercise.rds")
  
  # bar plots
  # healthy fibre
  sig_healthy_fibre_annot <- daa_healthy_fibre_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4 , !is.na(pathway_name))
  
  daa_healthy_fibre_annot_clean <- sig_healthy_fibre_annot |>
    mutate(
      pathway_class = pathway_class |>
        str_split(";") |>
        lapply(function(x) {
          x <- trimws(x)
          x <- x[-1]
          x <- str_replace(x, "^[0-9]+\\s*", "")
          paste(x, collapse = "; ")
        }) |>
        unlist()
    )
  
  daa_healthy_fibre_annot_clean$pathway_name <- gsub("\\s*\\[EC:[^]]+\\]$", "",daa_healthy_fibre_annot_clean$pathway_name)
  daa_healthy_fibre_annot_clean$pathway_class <-
    daa_healthy_fibre_annot_clean$pathway_class |>
    (\(x) gsub("^Unclassified;\\s*", "", x))() |>
    (\(x) gsub("^Not Included in Pathway or Brite;\\s*", "", x))()
  
  
  bar_health_fib <- ggplot(daa_healthy_fibre_annot_clean, aes(x = reorder(pathway_class, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#CC3C82", "#699CCC"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Healthy",
         x = NULL,
         y = "log2FC (fibre contrast)",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 13))
  
   #healthy exercise
  sig_healthy_exercise_annot <- daa_healthy_exercise_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4 , !is.na(pathway_name))
  
  daa_healthy_exercise_annot_clean <- sig_healthy_exercise_annot |>
    mutate(
      pathway_class = pathway_class |>
        str_split(";") |>
        lapply(function(x) {
          x <- trimws(x)
          x <- x[-1]
          x <- str_replace(x, "^[0-9]+\\s*", "")
          paste(x, collapse = "; ")
        }) |>
        unlist()
    )
  
  daa_healthy_exercise_annot_clean$pathway_name <- gsub("\\s*\\[EC:[^]]+\\]$", "",daa_healthy_exercise_annot_clean$pathway_name)
  daa_healthy_exercise_annot_clean$pathway_class <-
    daa_healthy_exercise_annot_clean$pathway_class |>
    (\(x) gsub("^Unclassified;\\s*", "", x))() |>
    (\(x) gsub("^Not Included in Pathway or Brite;\\s*", "", x))()
  
  
  bar_health_exercise <- ggplot(daa_healthy_exercise_annot_clean, aes(x = reorder(pathway_class, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#699CCC", "#CC3C82"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Healthy",
         x = NULL,
         y = "log2FC (exercise contrast)",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 13))
  
  # abnormal fibre
  sig_abnormal_fibre_annot <- daa_abnormal_fibre_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4 , !is.na(pathway_name))
  
  daa_abnormal_fibre_annot_clean <- sig_abnormal_fibre_annot |>
    mutate(
      pathway_class = pathway_class |>
        str_split(";") |>
        lapply(function(x) {
          x <- trimws(x)
          x <- x[-1]
          x <- str_replace(x, "^[0-9]+\\s*", "")
          paste(x, collapse = "; ")
        }) |>
        unlist()
    )
  
  daa_abnormal_fibre_annot_clean$pathway_name <- gsub("\\s*\\[EC:[^]]+\\]$", "",daa_abnormal_fibre_annot_clean$pathway_name)
  daa_abnormal_fibre_annot_clean$pathway_class <-
    daa_abnormal_fibre_annot_clean$pathway_class |>
    (\(x) gsub("^Unclassified;\\s*", "", x))() |>
    (\(x) gsub("^Not Included in Pathway or Brite;\\s*", "", x))()
  
  
  bar_abnormal_fib <- ggplot(daa_abnormal_fibre_annot_clean, aes(x = reorder(pathway_class, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#CC3C82", "#699CCC"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Abnormal",
         x = NULL,
         y = "log2FC (fibre contrast)",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 13))
  
  #abnormal exercise
  sig_abnormal_exercise_annot <- daa_abnormal_exercise_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4 , !is.na(pathway_name))
  
  daa_abnormal_exercise_annot_clean <- sig_abnormal_exercise_annot |>
    mutate(
      pathway_class = pathway_class |>
        str_split(";") |>
        lapply(function(x) {
          x <- trimws(x)
          x <- x[-1]
          x <- str_replace(x, "^[0-9]+\\s*", "")
          paste(x, collapse = "; ")
        }) |>
        unlist()
    )
  
  daa_abnormal_exercise_annot_clean$pathway_name <- gsub("\\s*\\[EC:[^]]+\\]$", "",daa_abnormal_exercise_annot_clean$pathway_name)
  daa_abnormal_exercise_annot_clean$pathway_class <-
    daa_abnormal_exercise_annot_clean$pathway_class |>
    (\(x) gsub("^Unclassified;\\s*", "", x))() |>
    (\(x) gsub("^Not Included in Pathway or Brite;\\s*", "", x))()
  
  
  bar_abnormal_exercise <- ggplot(daa_abnormal_exercise_annot_clean, aes(x = reorder(pathway_class, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#CC3C82", "#699CCC"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Abnormal",
         x = NULL,
         y = "log2FC (exercise contrast)",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none",axis.text.y = element_text(size = 13))
  
  # Top row = exercise
  top_row <- bar_health_exercise + bar_abnormal_exercise
  
  # Bottom row = fibre
  bottom_row <- bar_health_fib + bar_abnormal_fib
  
  # Combine
  bar_panel <- 
    top_row / 
    bottom_row + 
    plot_layout(heights = c(1, 1))  # shrink header row
  
  #Save the plot
  ggsave(plot = bar_panel, "results/aim4/05-pathway_class_Log2_Fold_Change_Plot.png", 
         width = 15, height = 9, units = "in", dpi = 300)
  
  # heatmap
  # healthy
  ko_healthy_relab <- ko_healthy |>
    apply(2, function(x) x/sum(x)) |> as.data.frame()
  
  sigs_healthy_fibre <- daa_healthy_fibre_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4) |>
    pull("feature") |> unique()
  
  sigs_healthy_exercise <- daa_healthy_exercise_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4) |>
    pull("feature") |> unique()
  
  stats_healthy_fibre <- ko_healthy_relab |> t() |> as.data.frame() |>
    select(all_of(sigs_healthy_fibre)) |>
    cor(method = "spearman")
  
  stats_healthy_exercise <- ko_healthy_relab |> t() |> as.data.frame() |>
    select(all_of(sigs_healthy_exercise)) |>
    cor(method = "spearman")
  
  pheatmap(stats_healthy_fibre,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",main = "Healthy CV, Fibre Contrast",
           filename = "results/aim4/06-healthy_fibre_heatmap.png")
  
  pheatmap(stats_healthy_exercise,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete", main = "Healthy CV, Exercise Contrast",
           filename = "results/aim4/07-healthy_exercise_heatmap.png")
  
  #abnormal
  ko_abnormal_relab <- ko_abnormal |>
    apply(2, function(x) x/sum(x)) |> as.data.frame()
  
  sigs_abnormal_fibre <- daa_abnormal_fibre_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4) |>
    pull("feature") |> unique()
  
  sigs_abnormal_exercise <- daa_abnormal_exercise_annot |>
    filter(p_adjust <= 0.01, abs(log2FoldChange) > 4) |>
    pull("feature") |> unique()
  
  stats_abnormal_fibre <- ko_abnormal_relab |> t() |> as.data.frame() |>
    select(all_of(sigs_abnormal_fibre)) |>
    cor(method = "spearman")
  
  stats_abnormal_exercise <- ko_abnormal_relab |> t() |> as.data.frame() |>
    select(all_of(sigs_abnormal_exercise)) |>
    cor(method = "spearman")
  
  pheatmap(stats_abnormal_fibre,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",main = "Abnormal CV, Fibre Contrast",
           filename = "results/aim4/08-abnormal_fibre_heatmap.png"
  )
  
  pheatmap(stats_abnormal_exercise,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",main = "Abnormal CV, Exercise Contrast",
           filename = "results/aim4/09-abnormal_exercise_heatmap.png"
  )
  
  ### pathwayPCA
  pathway_pca(abundance = ko_healthy, 
              metadata = meta_healthy,
              group = "fibre_group")
  
  ggsave("results/aim4/10-healthy_fibre_pca.png")
  
  pathway_pca(abundance = ko_healthy, 
              metadata = meta_healthy,
              group = "exercise_group",
              colors = c("#638475", "#f6d0b1"))
  
  ggsave("results/aim4/11-healthy_exercise_pca.png")
  
  pathway_pca(abundance = ko_abnormal, 
              metadata = meta_abnormal,
              group = "fibre_group")
  
  ggsave("results/aim4/12-abnormal_fibre_pca.png")
  
  pathway_pca(abundance = ko_abnormal, 
              metadata = meta_abnormal,
              group = "exercise_group",
              colors = c("#638475", "#f6d0b1"))
  
  ggsave("results/aim4/13-abnormal_exercise_pca.png")
  
  # save significant pathways
  write_tsv(daa_healthy_exercise_annot_clean, "results/aim4/14-sig_pathways_healthy_exercise.tsv")
  write_tsv(daa_healthy_fibre_annot_clean, "results/aim4/15-sig_pathways_healthy_fibre.tsv")
  write_tsv(daa_abnormal_exercise_annot_clean, "results/aim4/16-sig_pathways_abnormal_exercise.tsv")
  write_tsv(daa_abnormal_fibre_annot_clean, "results/aim4/17-sig_pathways_abnormal_fibre.tsv")
}

main()
