# This script performs differential abundance analysis for aim 3.
# Date: November 25th 2025
# Usage: Rscript R/06-aim3_DESeq.r

# Load the libraries
library(tidyverse)
library(tidymodels)
library(phyloseq)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(2025)

#read data
main <- function(){
  # Read in data
  phyloseq_obj <- readRDS("data/data_raw/phyloseq_object.rds")
  
  #rewrite metadata
  meta <- sample_data(phyloseq_obj)
  meta <- data.frame(meta)
  
  # Create the comparison columns
  meta <- meta |>
    mutate(fibre_group = ifelse(fiber >= 20, "high", "low"),
           exercise_group = ifelse(MET_mins_per_week > 1000, "high", "low"))
  
  meta$fibre_group <- as.factor(meta$fibre_group)
  meta$exercise_group <- as.factor(meta$exercise_group)
  
  meta$fibre_group <- relevel(meta$fibre_group, ref = "low")
  meta$exercise_group <- relevel(meta$exercise_group, ref = "low")
  
  meta$Cardiometabolic_status <- as.factor(meta$Cardiometabolic_status)
  meta$Cardiometabolic_status <- relevel(meta$Cardiometabolic_status, ref = "Healthy")
  
  sample_data(phyloseq_obj) <- meta
  
  #Add 1 to all counts because DESeq does not work well with zeros
  phyloseq_object_plus1 <- transform_sample_counts(phyloseq_obj, function(x) x+1)
  
  #separate
  phylo_healthy <- subset_samples(phyloseq_object_plus1, Cardiometabolic_status == "Healthy")
  phylo_abnormal <- subset_samples(phyloseq_object_plus1, Cardiometabolic_status == "Abnormal")
  
  
  ##### DESeq #####
  
  ############# healthy
  #Convert to DESeq dataset
  deseq_obj_healthy <- phyloseq_to_deseq2(
    phylo_healthy,
    ~ fibre_group + exercise_group + adiponectin + BMI + Calorie_intake + city + stool_consistency
  )
  
  #Run the DESeq2 model
  deseq_healthy <- DESeq(deseq_obj_healthy)
  
  #Extract results for fibre and exercise
  res_healthy_fibre <- results(deseq_healthy, contrast = c("fibre_group","high","low"))
  df_healthy_fibre <- as.data.frame(res_healthy_fibre)
  
  res_healthy_exercise <- results(deseq_healthy, contrast = c("exercise_group","high","low"))
  df_healthy_exercise <- as.data.frame(res_healthy_exercise)
  
  #Extract genus labels instead of just ASV IDs
  taxdf_healthy <- as.data.frame(tax_table(phylo_healthy))
  
  #append
  df_healthy_exercise_w_taxa <- cbind(df_healthy_exercise, taxdf_healthy)
  df_healthy_fibre_w_taxa <- cbind(df_healthy_fibre, taxdf_healthy)
  
  write_tsv(df_healthy_exercise_w_taxa, "results/aim3/DESeq2/01-healthy_exercise_res.tsv")
  write_tsv(df_healthy_fibre_w_taxa, "results/aim3/DESeq2/02-healthy_fibre_res.tsv")
  
  ############# abnormal
  #Convert to DESeq dataset
  deseq_obj_abnormal <- phyloseq_to_deseq2(
    phylo_abnormal,
    ~ fibre_group + exercise_group + adiponectin + BMI + Calorie_intake + city + stool_consistency
  )
  
  #Run the DESeq2 model
  deseq_abnormal <- DESeq(deseq_obj_abnormal)
  
  #Extract results for fibre and exercise
  res_abnormal_fibre <- results(deseq_abnormal, contrast = c("fibre_group","high","low"))
  df_abnormal_fibre <- as.data.frame(res_abnormal_fibre)
  
  res_abnormal_exercise <- results(deseq_abnormal, contrast = c("exercise_group","high","low"))
  df_abnormal_exercise <- as.data.frame(res_abnormal_exercise)
  
  #Extract genus labels instead of just ASV IDs
  taxdf_abnormal <- as.data.frame(tax_table(phylo_abnormal))
  
  #append
  df_abnormal_exercise_w_taxa <- cbind(df_abnormal_exercise, taxdf_abnormal)
  df_abnormal_fibre_w_taxa <- cbind(df_abnormal_fibre, taxdf_abnormal)
  
  write_tsv(df_abnormal_exercise_w_taxa, "results/aim3/DESeq2/03-abnormal_exercise_res.tsv")
  write_tsv(df_abnormal_fibre_w_taxa, "results/aim3/DESeq2/04-abnormal_fibre_res.tsv")
  
  ##### Plots #####                                             
  # make labels
  df_abnormal_exercise_w_taxa <- df_abnormal_exercise_w_taxa |>
    rowwise() |>
    mutate(label = make_label(c_across(Domain:Species))) |>
    ungroup()
  
  df_abnormal_fibre_w_taxa <- df_abnormal_fibre_w_taxa |>
    rowwise() |>
    mutate(label = make_label(c_across(Domain:Species))) |>
    ungroup()
  
  df_healthy_fibre_w_taxa <- df_healthy_fibre_w_taxa |>
    rowwise() |>
    mutate(label = make_label(c_across(Domain:Species))) |>
    ungroup()
  
  df_healthy_exercise_w_taxa <- df_healthy_exercise_w_taxa |>
    rowwise() |>
    mutate(label = make_label(c_across(Domain:Species))) |>
    ungroup()
  
  
  #Generating the volcano plots
  ab_ex_volcano <- EnhancedVolcano(as.data.frame(df_abnormal_exercise_w_taxa),
                  lab = df_abnormal_exercise_w_taxa$label, 
                  x = "log2FoldChange",
                  y = "padj",
                  pCutoff = 0.05,
                  FCcutoff = 1.5,
                  subtitle = NULL,
                  title = "Exercise Contrast",
                  legendPosition = "none",
                  caption = NULL,
                  axisLabSize = 10,
                  titleLabSize = 10,
                  labSize = 3,
                  pointSize = 2) + theme_classic() + theme(legend.position = "none")
  
  ab_fibre_volcano <- EnhancedVolcano(as.data.frame(df_abnormal_fibre_w_taxa),
                                   lab = df_abnormal_fibre_w_taxa$label, 
                                   x = "log2FoldChange",
                                   y = "padj",
                                   pCutoff = 0.05,
                                   FCcutoff = 1.5,
                                   subtitle = NULL,
                                   title = "Fibre Contrast",
                                   legendPosition = "none",
                                   caption = NULL,
                                   axisLabSize = 10,
                                   titleLabSize = 10,
                                   labSize = 3,
                                   pointSize = 2) + theme_classic() + theme(legend.position = "none")
  
  healthy_ex_volcano <- EnhancedVolcano(as.data.frame(df_healthy_exercise_w_taxa),
                                   lab = df_healthy_exercise_w_taxa$label, 
                                   x = "log2FoldChange",
                                   y = "padj",
                                   pCutoff = 0.05,
                                   FCcutoff = 1.5,
                                   subtitle = NULL,
                                   title = "Exercise Contrast",
                                   legendPosition = "none",
                                   caption = NULL,
                                   axisLabSize = 10,
                                   titleLabSize = 10,
                                   labSize = 3,
                                   pointSize = 2) + theme_classic() + theme(legend.position = "none")
  
  healthy_fibre_volcano <- EnhancedVolcano(as.data.frame(df_healthy_fibre_w_taxa),
                                   lab = df_healthy_fibre_w_taxa$label, 
                                   x = "log2FoldChange",
                                   y = "padj",
                                   pCutoff = 0.05,
                                   FCcutoff = 1.5,
                                   subtitle = NULL,
                                   title = "Fibre Contrast",
                                   legendPosition = "none",
                                   caption = NULL,
                                   axisLabSize = 10,
                                   titleLabSize = 10,
                                   labSize = 3,
                                   pointSize = 2) + theme_classic() + theme(legend.position = "none")
  
  title_healthy  <- ggplot() + 
    theme_void() + 
    ggtitle("Healthy CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  title_abnormal <- ggplot() + 
    theme_void() + 
    ggtitle("Abnormal CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Row 1: Exercise plots
  top_row <- healthy_ex_volcano + ab_ex_volcano
  
  # Row 2: Fibre plots
  bottom_row <- healthy_fibre_volcano + ab_fibre_volcano
  
  # Combine everything with column headers
  final_plot <- 
    (title_healthy + title_abnormal) /     # column titles
    (top_row) /                            # exercise row
    (bottom_row) +                         # fibre row
    plot_layout(heights = c(0.1, 1, 1))    # make titles smaller

  #Save the plots
  ggsave(plot = final_plot, "results/aim3/DESeq2/05-DESeq_Volcano_Plot.png", 
         width = 11, height = 12, units = "in", dpi = 300)
  
  #Generating bar plot
  df_healthy_fibre_w_taxa_sig <- df_healthy_fibre_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 2 | log2FoldChange <= -2))|>
    mutate(
      Family = gsub("^f__", "", Family),  
      Family = ifelse(grepl("Incertae_Sedis", Family, ignore.case = TRUE),
                      NA, Family)) |>
    filter(!is.na(Family))
  
  df_healthy_exercise_w_taxa_sig <- df_healthy_exercise_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 2 | log2FoldChange <= -2))|>
    mutate(
      Family = gsub("^f__", "", Family),  
      Phylum = gsub("^p__", "", Phylum),
      Family = ifelse(grepl("Incertae_Sedis", Family, ignore.case = TRUE),
                      NA, Family)) |>
    filter(!is.na(Family))
  
  df_abnormal_fibre_w_taxa_sig <- df_abnormal_fibre_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 2 | log2FoldChange <= -2))|>
    mutate(
      Family = gsub("^f__", "", Family),  
      Family = ifelse(grepl("Incertae_Sedis", Family, ignore.case = TRUE),
                      NA, Family)) |>
    filter(!is.na(Family))
  
  df_abnormal_exercise_w_taxa_sig <- df_abnormal_exercise_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 2 | log2FoldChange <= -2))|>
    mutate(
      Family = gsub("^f__", "", Family),  
      Family = ifelse(grepl("Incertae_Sedis", Family, ignore.case = TRUE),
                      NA, Family)) |>
    filter(!is.na(Family))
  
  bar_health_fib <- ggplot(df_healthy_fibre_w_taxa_sig, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Fibre Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_health_ex <- ggplot(df_healthy_exercise_w_taxa_sig, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Exercise Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_abnormal_fibre <- ggplot(df_abnormal_fibre_w_taxa_sig, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Fibre Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_abnormal_exercise <- ggplot(df_abnormal_exercise_w_taxa_sig, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Exercise Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  title_healthy  <- ggplot() + 
    theme_void() + 
    ggtitle("Healthy CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  title_abnormal <- ggplot() + 
    theme_void() + 
    ggtitle("Abnormal CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Top row = exercise
  top_row <- bar_health_ex + bar_abnormal_exercise
  
  # Bottom row = fibre
  bottom_row <- bar_health_fib + bar_abnormal_fibre
  
  # Combine
  bar_panel <- 
    (title_healthy + title_abnormal) /   # column headers
    top_row / 
    bottom_row +
    plot_layout(heights = c(0.1, 1, 1))  # shrink header row
  
  #Save the plot
  ggsave(plot = bar_panel, "results/aim3/DESeq2/06-DESeq_Log2_Fold_Change_Plot_family.png", 
         width = 10, height = 9, units = "in", dpi = 300)
  
  ### try with Phylum
  df_healthy_fibre_w_taxa_sig <- df_healthy_fibre_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))|>
    mutate(Phylum = gsub("^p__", "", Phylum))
  
  df_healthy_exercise_w_taxa_sig <- df_healthy_exercise_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))|>
    mutate(Phylum = gsub("^p__", "", Phylum))
  
  df_abnormal_fibre_w_taxa_sig <- df_abnormal_fibre_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))|>
    mutate(Phylum = gsub("^p__", "", Phylum))
  
  df_abnormal_exercise_w_taxa_sig <- df_abnormal_exercise_w_taxa |>
    filter(padj <= 0.05, (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))|>
    mutate(Phylum = gsub("^p__", "", Phylum))
  
  bar_health_fib <- ggplot(df_healthy_fibre_w_taxa_sig, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Fibre Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_health_ex <- ggplot(df_healthy_exercise_w_taxa_sig, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Exercise Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_abnormal_fibre <- ggplot(df_abnormal_fibre_w_taxa_sig, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Fibre Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  bar_abnormal_exercise <- ggplot(df_abnormal_exercise_w_taxa_sig, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # horizontal bars
    scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), labels = c("Under-represented","Over-represented")) +
    labs(title = "Exercise Contrast",
         x = NULL,
         y = "log2 Fold Change",
         fill = "Direction") +
    theme_classic()+
    theme(legend.position = "none")
  
  title_healthy  <- ggplot() + 
    theme_void() + 
    ggtitle("Healthy CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  title_abnormal <- ggplot() + 
    theme_void() + 
    ggtitle("Abnormal CV") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Top row = exercise
  top_row <- bar_health_ex + bar_abnormal_exercise
  
  # Bottom row = fibre
  bottom_row <- bar_health_fib + bar_abnormal_fibre
  
  # Combine
  bar_panel <- 
    (title_healthy + title_abnormal) /   # column headers
    top_row / 
    bottom_row +
    plot_layout(heights = c(0.1, 1, 1))  # shrink header row
  
  #Save the plot
  ggsave(plot = bar_panel, "results/aim3/DESeq2/07-DESeq_Log2_Fold_Change_Plot_phylum.png", 
         width = 7, height = 5, units = "in", dpi = 300)
  
  ############################################## interactions
  #Convert to DESeq dataset
  deseq_obj_fibre <- phyloseq_to_deseq2(
    phyloseq_object_plus1,
    ~ Cardiometabolic_status + fibre_group + Cardiometabolic_status:fibre_group + BMI + Calorie_intake + city + stool_consistency
  )
  
  deseq_obj_exercise <- phyloseq_to_deseq2(
    phyloseq_object_plus1,
    ~ Cardiometabolic_status + exercise_group + Cardiometabolic_status:exercise_group + BMI + Calorie_intake + city + stool_consistency
  )
  
  #Run the DESeq2 model
  deseq_fibre <- DESeq(deseq_obj_fibre)
  deseq_exercise <- DESeq(deseq_obj_exercise)
  
  #Extract results for fibre and exercise
  res_cv_fibre <- results(deseq_fibre, name = "Cardiometabolic_statusAbnormal.fibre_grouphigh")
  df_cv_fibre <- as.data.frame(res_cv_fibre)
  
  res_cv_exercise <- results(deseq_exercise, name = "Cardiometabolic_statusAbnormal.exercise_grouphigh")
  df_cv_exercise <- as.data.frame(res_cv_exercise)
  

  #Extract genus labels instead of just ASV IDs
  taxdf <- as.data.frame(tax_table(phyloseq_object_plus1))
  
  #append
  df_cv_exercise_w_taxa <- cbind(df_cv_exercise, taxdf)
  df_cv_fibre_w_taxa <- cbind(df_cv_fibre, taxdf)
  
  write_tsv(df_cv_exercise_w_taxa, "results/aim3/DESeq2/08-exercise_cv_interaction_res.tsv")
  write_tsv(df_cv_fibre_w_taxa, "results/aim3/DESeq2/09-cv_fibre_cv_interaction_res.tsv")
 }

# helper functions
make_label <- function(x) {
  
  # 1. If Genus exists AND is not Incertae_Sedis
  if (!is.na(x["Genus"]) && x["Genus"] != "" && !grepl("Incertae_Sedis", x["Genus"], ignore.case = TRUE)) {
    
    # If Species exists AND is not Incertae_Sedis → return Genus + Species
    if (!is.na(x["Species"]) && x["Species"] != "" &&
        !grepl("Incertae_Sedis", x["Species"], ignore.case = TRUE)) {
      
      return(paste(x["Genus"], x["Species"]))
    }
    
    # Otherwise return Genus only
    return(x["Genus"])
  }
  
  # 2. Otherwise fallback: find deepest non–Incertae_Sedis rank
  cleaned <- x[!is.na(x) & x != "" & !grepl("Incertae_Sedis", x, ignore.case = TRUE)]
  
  if (length(cleaned) > 0) {
    deepest <- tail(cleaned, 1)
    return(as.character(deepest))
  }
  
  # 3. If all ranks missing or Incertae_Sedis
  return(NA_character_)
}

main()
