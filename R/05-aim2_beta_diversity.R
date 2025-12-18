# This script performs beta diversity under aim 2. 
# Date: November 16, 2025
# Usage: Rscript R/05-aim2_beta_diversity.R

# Load in libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(picante)
library(ggplot2)
library(cowplot)
library(viridis)
library(tidyr)
library(RColorBrewer)
library(patchwork)
# Set the seed for reproducibility 
set.seed(2025)

main <- function(){
# Read in data
  phyloseq_obj <- readRDS("data/data_processed/phyloseq_obj_rarefied.rds")
  
  #rewrite metadata
  meta <- sample_data(phyloseq_obj)
  meta <- data.frame(meta)

  # Create the lifestyle_group column
  meta <- meta |>
    mutate(
      lifestyle_group = case_when(
        fiber >= 20 & MET_mins_per_week > 1000 ~ "adequate fibre;high exercise",
        fiber >= 20 & MET_mins_per_week <= 1000  ~ "adequate fibre;low exercise",
        fiber < 20  & MET_mins_per_week > 1000 ~ "inadequate fibre;high exercise",
        TRUE ~ "inadequate fibre;low exercise"
      )
    )
  
  meta$group <- paste(meta$lifestyle_group, meta$Cardiometabolic_status, sep = "_")
  
  sample_data(phyloseq_obj) <- meta
  
  # subset
  phylo_healthy <- subset_samples(phyloseq_obj, Cardiometabolic_status == "Healthy")
  phylo_abnormal <- subset_samples(phyloseq_obj, Cardiometabolic_status == "Abnormal")
  
  #PCoA ordination and permanova
  dm_braycurtis_healthy <- vegdist(t(otu_table(phylo_healthy)), method="bray")
  dm_braycurtis_abnormal <- vegdist(t(otu_table(phylo_abnormal)), method="bray")
  
  saveRDS(dm_braycurtis_healthy, "results/aim2/beta_diversity/01-dm_braycurtis_healthy.rds")
  saveRDS(dm_braycurtis_abnormal, "results/aim2/beta_diversity/02-dm_braycurtis_abnormal.rds")
  
  permanova_healthy <- adonis2(dm_braycurtis_healthy ~ lifestyle_group, data=data.frame(sample_data(phylo_healthy)))
  permanova_abnormal <- adonis2(dm_braycurtis_abnormal ~ lifestyle_group, data=data.frame(sample_data(phylo_abnormal)))
  
  saveRDS(permanova_healthy, "results/aim2/beta_diversity/03-permanova_healthy.rds")
  saveRDS(permanova_abnormal, "results/aim2/beta_diversity/04-permanova_abnormal.rds")
  
  #PCoA
  ord_healthy <- ordinate(phylo_healthy, method = "PCoA", distance = dm_braycurtis_healthy)
  ord_abnormal <- ordinate(phylo_abnormal, method = "PCoA", distance = dm_braycurtis_abnormal)
  
  saveRDS(ord_healthy, "results/aim2/beta_diversity/05-ordination_healthy.rds")
  saveRDS(ord_abnormal, "results/aim2/beta_diversity/06-ordination_abnormal.rds")
  
  # plot
  pval_healthy <- permanova_healthy$`Pr(>F)`[1]
  pval_abnormal <- permanova_abnormal$`Pr(>F)`[1]
  
  var_explained_healthy <- ord_healthy$values$Relative_eig * 100
  PC1_pct_healthy <- round(var_explained_healthy[1], 1)
  PC2_pct_healthy <- round(var_explained_healthy[2], 1)
  
  bray_plot_healthy <- plot_ordination(phylo_healthy, ord_healthy, color = "lifestyle_group") +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2) + 
    labs(x = paste0("PCoA1 (", PC1_pct_healthy, "%)"),
         y = paste0("PCoA2 (", PC2_pct_healthy, "%)"),
         title = "Healthy",
         subtitle = paste0("PERMANOVA p = ", pval_healthy),
         color = "Lifestyle \ngroup") +
    scale_color_brewer(
      palette = "Set2",
      labels = c(
        "adequate fibre;high exercise"    = "High fibre\nHigh exercise",
        "adequate fibre;low exercise"     = "High fibre\nLow exercise",
        "inadequate fibre;high exercise"  = "Low fibre\nHigh exercise",
        "inadequate fibre;low exercise"   = "Low fibre\nLow exercise"
      )) +
    theme_classic()+
    theme(legend.position = "none",  axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          strip.text = element_text(size = 13, face = "bold"))
  
  var_explained_abnormal <- ord_abnormal$values$Relative_eig * 100
  PC1_pct_abnormal <- round(var_explained_abnormal[1], 1)
  PC2_pct_abnormal <- round(var_explained_abnormal[2], 1)
  
  bray_plot_abnormal <- plot_ordination(phylo_abnormal, ord_abnormal, color = "lifestyle_group") +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2) + 
    labs(x = paste0("PCoA1 (", PC1_pct_abnormal, "%)"),
         y = paste0("PCoA2 (", PC2_pct_abnormal, "%)"),
         title = "Abnormal",
         subtitle = paste0("PERMANOVA p = ", pval_abnormal),
         color = "Lifestyle \ngroup") +
    scale_color_brewer(
      palette = "Set2",
      labels = c(
        "adequate fibre;high exercise"    = "High fibre\nHigh exercise",
        "adequate fibre;low exercise"     = "High fibre\nLow exercise",
        "inadequate fibre;high exercise"  = "Low fibre\nHigh exercise",
        "inadequate fibre;low exercise"   = "Low fibre\nLow exercise"
      )) +
    theme_classic()+
    theme(legend.position = "none",  axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          axis.title.x = element_text(size = 13),
          strip.text = element_text(size = 13, face = "bold"))
  
  # combine plots
  final_plot <- (bray_plot_healthy + bray_plot_abnormal) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 14))
  
  ggsave(plot = final_plot, "results/aim2/beta_diversity/07-pcoa_plot.png",
         width = 10, height = 6)
}

main()
