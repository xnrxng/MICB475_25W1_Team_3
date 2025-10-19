# This script performs preliminary exploratory data analysis on the data.
# Date: October 19th 2025
# Usage: Rscript R/01-prelim_eda.R

library(tidyverse)
library(utils)
library(data.table)
library(ggplot2)
library(cowplot)
library(vegan)

main <- function(){
  # read in metadata
  meta <- fread("data/data_raw/colombia_metadata.txt")
  
  # make histogram of calorie intake distribution colored by sex and BMI category. add median lines
  sex_median_values_cal <- meta |>
    group_by(sex) |>
    summarize(median_calorie = median(Calorie_intake, na.rm = TRUE))
  
  bmi_median_values_ca <- meta |>
    group_by(BMI_class) |>
    summarize(median_calorie = median(Calorie_intake, na.rm = TRUE))
  
  
  p_cal_sex <- ggplot(meta, aes(x = Calorie_intake, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily calorie intake (kcal)",
      y = "Number of people",
      title = "Calorie Intake Distribution per Sex",
      color = "Biological Sex",
      fill = "Biological Sex"
    ) +
    theme_classic() +
    scale_fill_manual(
      values = c("male" = "#df8e5f", "female" = "#67dce5")
    ) +
    geom_vline(
      data = sex_median_values_cal,
      aes(xintercept = median_calorie, color = sex),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("male" = "#b35e27", "female" = "#1fa2b3")
    )+
    theme(legend.position = "none")
  
  p_cal_bmi <- ggplot(meta, aes(x = Calorie_intake, fill = BMI_class, color = BMI_class)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Daily calorie intake (kcal)",
      y = "Number of people",
      title = "Calorie Intake Distribution per BMI Category",
      color = "BMI Category",
      fill = "BMI Category"
    ) +
    scale_color_manual(
      values = c("Overweight" = "#0D6E1F", "Obese" = "#C20010", "Lean" = "#502682")
    ) +
    scale_fill_manual(
      values = c("Overweight" = "#B4DC7F", "Obese" = "#FFA0AC", "Lean" = "#C9B1E7")
    ) +
    geom_vline(
      data = bmi_median_values_ca,
      aes(xintercept = median_calorie, color = BMI_class),
      linetype = "dashed",
      linewidth = 1
    ) +
    theme_classic()+
    theme(legend.position = "none")
  
  # make histogram of calorie intake distribution colored by sex and BMI category. add median lines
  sex_median_values <- meta |>
    group_by(sex) |>
    summarize(median_age = median(age_years, na.rm = TRUE))
  
  bmi_median_values <- meta |>
    group_by(BMI_class) |>
    summarize(median_age = median(age_years, na.rm = TRUE))
  
  p_age_sex <- ggplot(meta, aes(x = age_years, fill = sex, color = sex)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Age (years)",
      y = "Number of people",
      title = "Age Distribution per Sex",
      color = "Biological Sex",
      fill = "Biological Sex"
    ) +
    theme_classic() +
    scale_fill_manual(
      values = c("male" = "#df8e5f", "female" = "#67dce5")
    ) +
    geom_vline(
      data = sex_median_values,
      aes(xintercept = median_age, color = sex),
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("male" = "#b35e27", "female" = "#1fa2b3")
    )
  
  p_age_bmi <- ggplot(meta, aes(x = age_years, fill = BMI_class, color = BMI_class)) +
    geom_histogram(alpha = 0.6) +
    labs(
      x = "Age (years)",
      y = "Number of people",
      title = "Calorie Intake Distribution per BMI Category",
      color = "BMI Category",
      fill = "BMI Category"
    ) +
    theme_classic() +
    scale_color_manual(
      values = c("Overweight" = "#0D6E1F", "Obese" = "#C20010", "Lean" = "#502682")
    ) +
    scale_fill_manual(
      values = c("Overweight" = "#B4DC7F", "Obese" = "#FFA0AC", "Lean" = "#C9B1E7")
    ) +
    geom_vline(
      data = bmi_median_values,
      aes(xintercept = median_age, color = BMI_class),
      linetype = "dashed",
      linewidth = 1
    )
  
  plotlist_hist <- list(p_cal_sex, p_age_sex, p_cal_bmi, p_age_bmi)
  
  # combine 4 plots into one and save
  comb_hists <- cowplot::plot_grid(plotlist = plotlist_hist, ncol = 2)
  cowplot::save_plot(plot = comb_hists, "results/01-bmi_and_age_hists.png", base_height = 8, base_width = 10)
  
  # remove outliers which is ppl with over 3k calories and under 1.2k calories
  meta_filt <- meta[(meta$Calorie_intake < 3000) & (meta$Calorie_intake > 1200), , drop = FALSE]
  
  # see number of healthy vs abnormal status
  counts_cardiometabolic_status <- meta_filt |> group_by(Cardiometabolic_status) |> count()
  saveRDS(counts_cardiometabolic_status, "results/02-counts_cardiometabolic_status.rds")
  
}

# helper functions

main()