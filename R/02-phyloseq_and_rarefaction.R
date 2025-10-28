# This script converts the files into a phyloseq object and rarefies it.
# Date: October 20th 2025
# Usage: Rscript R/02-phyloseq_and_rarefaction.R

library(tidyverse)
library(readr)
library(utils)
library(data.table)
library(ggplot2)
library(cowplot)
library(vegan)
library(Matrix)
library(stats)
library(ape)
library(phyloseq)
set.seed(2025)

main <- function(){
  meta <- read.delim("data/data_raw/colombia_metadata.txt", sep = "\t")
  otu <- read.delim(file = "data/data_processed/feature-table.txt", sep = "\t", skip = 1, row.names = 1)
  taxa <- read.delim(file = "data/data_processed/taxonomy.tsv", sep = "\t")
  tree <- read.tree(file = "data/data_processed/tree.nwk")
  
  ### converting into df, changing rownames to be the sample names, getting rid of sample_name
  meta_df <- as.data.frame(meta[,-1])
  rownames(meta_df) <- meta$X.SampleID
  sample <- sample_data(meta_df)
  
  ### converting it into matrix 
  otu_mat <- as.matrix(otu)
  otu_tb <- otu_table(otu_mat, taxa_are_rows = TRUE)
  
  ### getting rid of confidence column, separating taxon into 7 columns, converting into matrix
  taxa_mat <- taxa |> 
    dplyr::select(-Confidence) |>
    separate(col = Taxon, sep = "; ",
             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) |>
    as.matrix()
  
  ### getting rid of feature id as a column and making it the rownames
  taxa_mat <- taxa_mat[,-1]
  rownames(taxa_mat) <- taxa$`Feature.ID`
  taxa_tb <- tax_table(taxa_mat)
  
  ### converting it into the phyloseq object
  phyloseq_obj <- phyloseq(otu_tb, sample, taxa_tb, tree)
  
  ### keeping bacteria, getting rid of chloroplasts and mitochondria
  phyloseq_obj <- subset_taxa(phyloseq_obj, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
  saveRDS(phyloseq_obj, "data/data_raw/phyloseq_object.rds")
  
  ### rarecurve plot and rarefying
  sample_ids <- colnames(otu_mat)
  colors <- rainbow(length(sample_ids))
  
  png("results/06-rarefaction_curve.png")
  
  rarecurve(
    t(as.data.frame(otu_table(phyloseq_obj))),
    step = 20,
    col = colors,
    label = FALSE,                  
    cex = 0.6
  )
  
  abline(v = 24406, col = "red", lty = 2, lwd = 2)
  
  dev.off()

  phyloseq_obj_rare <- rarefy_even_depth(phyloseq_obj, sample.size = 24406, rngseed = 2025, verbose = TRUE)
  saveRDS(phyloseq_obj_rare, "data/data_processed/phyloseq_obj_rarefied.rds")
  
  ### remake hist plots with fewer samples
  meta_rare <- sample_data(phyloseq_obj_rare)
  
  # make histogram of fiber intake distribution colored by sex and cardiovascular health category. add median lines
  sex_median_values_fiber <- meta_rare |>
    group_by(sex) |>
    summarize(median_fiber = median(fiber, na.rm = TRUE))
  
  cv_median_values_fiber <- meta_rare |>
    group_by(Cardiometabolic_status) |>
    summarize(median_fiber = median(fiber, na.rm = TRUE))
  
  
  p_fiber_sex <- ggplot(meta_rare, aes(x = fiber, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily fibre intake (g)",
      y = "Number of people",
      title = "Fibre Intake Distribution per Sex",
      color = "Biological Sex",
      fill = "Biological Sex"
    ) +
    theme_classic() +
    scale_fill_manual(
      values = c("male" = "#df8e5f", "female" = "#67dce5")
    ) +
    geom_vline(
      data = sex_median_values_fiber,
      aes(xintercept = median_fiber, color = sex),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("male" = "#b35e27", "female" = "#1fa2b3")
    )+
    theme(legend.position = "none")
  
  p_fiber_cv <- ggplot(meta_rare, aes(x = fiber, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily fibre intake (g)",
      y = "Number of people",
      title = "Fibre Intake Distribution per \nCardiometabolic Status",
      color = "Cardiometabolic \nstatus",
      fill = "Cardiometabolic \nstatus"
    ) +
    scale_color_manual(
      values = c("Healthy" = "#0D6E1F", "Abnormal" = "#C20010")
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#B4DC7F", "Abnormal" = "#FFA0AC")
    ) +
    geom_vline(
      data = cv_median_values_fiber,
      aes(xintercept = median_fiber, color = Cardiometabolic_status),
      linetype = "dashed",
      linewidth = 1
    ) +
    theme_classic()+
    theme(legend.position = "none")
  
  # make histogram of met distribution colored by sex and cv category. add median lines
  sex_median_values_met <- meta_rare |>
    group_by(sex) |>
    summarize(median_met = median(MET_mins_per_week, na.rm = TRUE))
  
  cv_median_values_met <- meta_rare |>
    group_by(Cardiometabolic_status) |>
    summarize(median_met = median(MET_mins_per_week, na.rm = TRUE))
  
  p_met_sex <- ggplot(meta_rare, aes(x = MET_mins_per_week, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Metabolic Equivalent of Task (min/week)",
      y = "Number of people",
      title = "MET Distribution per Sex",
      color = "Biological Sex",
      fill = "Biological Sex"
    ) +
    theme_classic() +
    scale_fill_manual(
      values = c("male" = "#df8e5f", "female" = "#67dce5")
    ) +
    geom_vline(
      data = sex_median_values_met,
      aes(xintercept = median_met, color = sex),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("male" = "#b35e27", "female" = "#1fa2b3")
    )+
    theme(legend.position = "none")
  
  p_met_cv <- ggplot(meta_rare, aes(x = MET_mins_per_week, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Metabolic Equivalent of Task (min/week)",
      y = "Number of people",
      title = "MET Distribution per \nCardiometabolic status"
    ) +
    theme_classic() +
    scale_color_manual(
      values = c("Healthy" = "#0D6E1F", "Abnormal" = "#C20010")
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#B4DC7F", "Abnormal" = "#FFA0AC")
    ) +
    geom_vline(
      data = cv_median_values_met,
      aes(xintercept = median_met, color = Cardiometabolic_status),
      linetype = "dashed",
      linewidth = 1
    )+
    theme(legend.position = "none")
  
  # make histogram of adiponectin distribution colored by sex and cardiovascular health category. add median lines
  sex_median_values_adi <- meta_rare |>
    group_by(sex) |>
    summarize(median_adi = median(adiponectin, na.rm = TRUE))
  
  cv_median_values_adi <- meta_rare |>
    group_by(Cardiometabolic_status) |>
    summarize(median_adi = median(adiponectin, na.rm = TRUE))
  
  
  p_adi_sex <- ggplot(meta_rare, aes(x = adiponectin, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Adiponectin level (\U00B5g/mL)",
      y = "Number of people",
      title = "Adiponectin Distribution per Sex",
      color = "Biological Sex",
      fill = "Biological Sex"
    ) +
    theme_classic() +
    scale_fill_manual(
      values = c("male" = "#df8e5f", "female" = "#67dce5")
    ) +
    geom_vline(
      data = sex_median_values_adi,
      aes(xintercept = median_adi, color = sex),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("male" = "#b35e27", "female" = "#1fa2b3")
    )
  
  p_adi_cv <- ggplot(meta_rare, aes(x = adiponectin, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Adiponectin level (\U00B5g/mL)",
      y = "Number of people",
      title = "Adiponectin Distribution per \nCardiometabolic Status",
      color = "Cardiometabolic \nstatus",
      fill = "Cardiometabolic \nstatus"
    ) +
    scale_color_manual(
      values = c("Healthy" = "#0D6E1F", "Abnormal" = "#C20010")
    ) +
    scale_fill_manual(
      values = c("Healthy" = "#B4DC7F", "Abnormal" = "#FFA0AC")
    ) +
    geom_vline(
      data = cv_median_values_adi,
      aes(xintercept = median_adi, color = Cardiometabolic_status),
      linetype = "dashed",
      linewidth = 1
    ) +
    theme_classic()
  
  # combine 6 plots into one and save
  plotlist_hist <- list(p_fiber_sex, p_met_sex, p_adi_sex, p_fiber_cv, p_met_cv, p_adi_cv)
  comb_hists <- cowplot::plot_grid(plotlist = plotlist_hist, ncol = 3)
  cowplot::save_plot(plot = comb_hists, "results/07-fiber_met_adi_hists_after_rarefaction.png", base_height = 6, base_width = 12)
  
  # count number of patients that fall into each group
  meta_rare <- data.frame(meta_rare)
  meta_rare_mutated <- meta_rare |>
    mutate(fibre_group = ifelse(fiber >= 20, "high", "low"),
           exercise_group = ifelse(MET_mins_per_week > 1000, "high", "low"))
  
  meta_counts <- meta_rare_mutated |>
    group_by(Cardiometabolic_status, fibre_group, exercise_group) |>
    count()
  
  write_tsv(meta_counts, "results/08-cvstatus_fibre_exercise_counts.tsv")
}

#helper functions
main()