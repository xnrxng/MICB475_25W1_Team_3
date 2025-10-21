# This script performs preliminary exploratory data analysis on the data.
# Date: October 19th 2025
# Usage: Rscript R/01-prelim_eda.R

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
  # read in metadata
  meta <- fread("data/data_raw/colombia_metadata.txt")
  
  # remove outliers which is ppl with over 3k calories and under 1.2k calories
  #meta_filt <- meta[(meta$Calorie_intake < 3000) & (meta$Calorie_intake > 1200), , drop = FALSE]
  # remove rows with NA values
  meta_filt <- na.omit(meta)
  
  # make rownames sample id
  meta_filt <- as.data.frame(meta_filt)
  rownames(meta_filt) <- meta_filt$`#SampleID`
  meta_filt <- meta_filt |>
    dplyr::select(-`#SampleID`)
  
  # summary table of metadata
  summary_table <- meta_filt |>
    summarise(across(everything(), ~ {
      if (is.numeric(.x)) {
        paste0(min(.x, na.rm = TRUE), " – ", max(.x, na.rm = TRUE))
      } else {
        vals <- unique(.x)
        paste(sort(vals), collapse = ", ")
      }
    })) |>
    pivot_longer(cols = everything(),
                 names_to = "Variable",
                 values_to = "Range_or_Categories")
  write_tsv(summary_table, "results/0-metadata_summary.tsv")
  
  # get rid of country since it's only one category. also latitude since it is redundant with city
  meta_filt <- meta_filt |>
    select(-country, -latitude)
  
  # make histogram of fiber intake distribution colored by sex and cardiovascular health category. add median lines
  sex_median_values_fiber <- meta |>
    group_by(sex) |>
    summarize(median_fiber = median(fiber, na.rm = TRUE))
  
  cv_median_values_fiber <- meta |>
    group_by(Cardiometabolic_status) |>
    summarize(median_fiber = median(fiber, na.rm = TRUE))
  
  
  p_fiber_sex <- ggplot(meta, aes(x = fiber, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily fiber intake (g)",
      y = "Number of people",
      title = "Fiber Intake Distribution per Sex",
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
  
  p_fiber_cv <- ggplot(meta, aes(x = fiber, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily fiber intake (g)",
      y = "Number of people",
      title = "Fiber Intake Distribution per Cardiometabolic Status",
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
  sex_median_values_met <- meta |>
    group_by(sex) |>
    summarize(median_met = median(MET_mins_per_week, na.rm = TRUE))
  
  cv_median_values_met <- meta |>
    group_by(Cardiometabolic_status) |>
    summarize(median_met = median(MET_mins_per_week, na.rm = TRUE))
  
  p_met_sex <- ggplot(meta, aes(x = MET_mins_per_week, fill = sex, color = sex)) +
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
  
  p_met_cv <- ggplot(meta, aes(x = MET_mins_per_week, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Metabolic Equivalent of Task (min/week)",
      y = "Number of people",
      title = "MET Distribution per Cardiometabolic_status"
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
  sex_median_values_adi <- meta |>
    group_by(sex) |>
    summarize(median_adi = median(adiponectin, na.rm = TRUE))
  
  cv_median_values_adi <- meta |>
    group_by(Cardiometabolic_status) |>
    summarize(median_adi = median(adiponectin, na.rm = TRUE))
  
  
  p_adi_sex <- ggplot(meta, aes(x = adiponectin, fill = sex, color = sex)) +
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
  
  p_adi_cv <- ggplot(meta, aes(x = adiponectin, fill = Cardiometabolic_status, color = Cardiometabolic_status)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Adiponectin level (\U00B5g/mL)",
      y = "Number of people",
      title = "Adiponectin Distribution per Cardiometabolic Status",
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
  cowplot::save_plot(plot = comb_hists, "results/01-fiber_met_adi_hists.png", base_height = 6, base_width = 12)
  
  # see number of healthy vs abnormal status
  counts_cardiometabolic_status <- meta |> group_by(Cardiometabolic_status) |> count()
  write_tsv(counts_cardiometabolic_status, "results/02-counts_cardiometabolic_status.tsv")
  
  # read otu table 
  otu_table <- fread("data/data_processed/feature-table.txt")
  
  # bulk counts based on phylum
  taxonomy_tbl <- read_tsv("data/data_processed/taxonomy.tsv")
  taxonomy_mat <- taxonomy_tbl |> 
    dplyr::select(-Confidence) |>
    separate(col = Taxon, sep = "; ",
             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  
  merge_tbl <- merge(otu_table, taxonomy_mat, by.x = "#OTU ID", by.y = "Feature ID")
  
  otu_phylum <- merge_tbl |>
    group_by(Phylum) |>
    summarise(across(where(is.numeric), sum))
  
  otu_phylum_clean <- otu_phylum[, colnames(otu_phylum) != meta$`#SampleID`[!(meta$`#SampleID` %in% rownames(meta_filt))]]
  otu_phylum_clean <- na.omit(otu_phylum_clean)
  otu_phylum_clean <- as.data.frame(otu_phylum_clean)
  rownames(otu_phylum_clean) <- otu_phylum_clean$Phylum
  otu_phylum_clean <- otu_phylum_clean |>
    select(-Phylum)
  
  # redundancy analysis
  otu_phylum_clean <- otu_phylum_clean[, rownames(meta_filt)]
  otu_phylum_clean[] <- lapply(otu_phylum_clean, as.numeric)
  otu_phylum_clean_t <- t(otu_phylum_clean)
  
  otu_hel <- decostand(otu_phylum_clean_t, method = "hellinger")

  rda_model <- vegan::rda(otu_hel ~ . , data = meta_filt)
  
  saveRDS(rda_model, "results/03-prelim_rda_model.rds")
  
  anova_rda <- anova(rda_model, by = "term", permutations = 999)
  anova_rda <- rownames_to_column(anova_rda, var = "variable")
  
  write_tsv(anova_rda, "results/04-prelim_anova_rda.tsv")
  
  # prelim logistic regression of cardiometabolic status based on metadata only
  meta_filt$Cardiometabolic_status <- factor(meta_filt$Cardiometabolic_status,
                                             levels = c("Healthy", "Abnormal"))
  
  mdl <- glm(Cardiometabolic_status ~ ., data = meta_filt, family = binomial)
  mdl_summary <- summary(mdl)$coefficients
  mdl_summary <- as.data.frame(mdl_summary)
  mdl_summary$Variable <- rownames(mdl_summary)
  write_tsv(mdl_summary, "results/05-log_reg_summary.tsv")
}

# helper functions

main()
