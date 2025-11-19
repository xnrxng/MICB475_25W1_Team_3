# This script fits the 3 models from aim 1 onto the data.
# Date: November 3rd 2025
# Usage: Rscript R/03-aim1_predictors.R

library(tidyverse)
library(tidymodels)
library(readr)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(GGally)
library(glmnet)
library(glmbb)
library(broom)
library(leaps)
library(faraway)
library(mltools)
library(pROC)
library(caret)
library(car)
library(vegan)
library(ape)
library(phyloseq)
library(pheatmap)
library(ggordiplots)
set.seed(2025)

main <- function(){
  #read data
  phylo_rare <- readRDS("data/data_processed/phyloseq_obj_rarefied.rds")
  
  # aggregate
  ps_phylum <- tax_glom(phylo_rare, taxrank = "Phylum")
  
  ############ MLR
  # TSS transform
  ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
  
  # arcsine transform
  ps_phylum_asin <- transform_sample_counts(ps_phylum_rel, function(x) asin(sqrt(x)))
  
  # rename rownames to phylums and transpose
  phylum_asin <- otu_table(ps_phylum_asin)
  tax_table_data <- tax_table(ps_phylum_asin)
  rownames(phylum_asin) <- tax_table_data[, "Phylum"] |> as.vector()
  phylum_asin <- t(phylum_asin)
  phylum_asin <- as.data.frame(phylum_asin)
  phylum_asin$sample_id <- rownames(phylum_asin)
  
  # merge with metadata
  meta <- sample_data(ps_phylum_asin)
  meta <- data.frame(meta)
  meta$sample_id <- rownames(meta)
  meta <- meta |> select(sample_id, Cardiometabolic_status)
  merged_df <- merge(phylum_asin, meta, by = "sample_id")
  merged_df <- merged_df |> select(-sample_id)
  
  # convert CV to binary and run MLR
  merged_df <- merged_df |>
    mutate(Cardiometabolic_status = ifelse(Cardiometabolic_status == "Healthy", 0, 1))
  
  mlr_model <- lm(Cardiometabolic_status ~ ., data = merged_df)
  
  #save results
  saveRDS(mlr_model, "results/aim1/01-mlr_model.rds")
  summary_df <- as.data.frame(summary(mlr_model)$coefficients)
  summary_df$variable <- rownames(summary_df)
  write_tsv(summary_df, "results/aim1/02-mlr_coefficients.tsv")
  
  #############RDA
  # prepare counts for RDA
  phylum_tb <- otu_table(ps_phylum)
  tax_table_data <- tax_table(ps_phylum)
  rownames(phylum_tb) <- tax_table_data[, "Phylum"] |> as.vector()
  phylum_tb <- as.data.frame(phylum_tb)
  
  # remove country bc it's only 1 level, latitude and age range as well bc it's redundant
  meta <- sample_data(ps_phylum)
  meta <- data.frame(meta)
  meta <- meta |> select(-country, -latitude, -age_range)
  
  # make sure rows are in the same order
  phylum_tb <- t(phylum_tb)
  phylum_tb <- phylum_tb[rownames(meta), , drop = FALSE]
  
  # hellinger transform and run RDA
  otu_hel <- decostand(phylum_tb, method = "hellinger")
  rda_model <- vegan::rda(otu_hel ~ . , data = meta)
  
  saveRDS(rda_model, "results/aim1/03-rda_model.rds")
  
  # test significance
  anova_rda <- anova(rda_model, by = "term", permutations = 999)
  anova_rda <- rownames_to_column(anova_rda, var = "variable")
  
  write_tsv(anova_rda, "results/aim1/04-rda_coefficients.tsv")
  
  #refit significant variables with adonis2
  adonis <- adonis2(otu_hel ~ city + stool_consistency+adiponectin + Calorie_intake + BMI + fiber, data = meta, by = "terms")
  adonis <- rownames_to_column(adonis, var = "variable")
  write_tsv(adonis, "results/aim1/05-adonis_sig_rda.tsv")
  
  # ordination plots
  meta <- meta |>
    mutate(fibre_group = ifelse(fiber >= 20, "high", "low"),
           exercise_group = ifelse(MET_mins_per_week > 1000, "high", "low"),
           fibre_exercise_group = paste(fibre_group, exercise_group, sep = "_"))
  
  ggordi <- gg_ordiplot(
    ord = rda_model, 
    groups = meta$fibre_exercise_group, 
    kind = "se", 
    conf = 0.95,
    pt.size = 2, 
    plot = TRUE
  )
  
  p <- ggordi$plot +
    labs(
      title = "Phylum - RDA",
      subtitle = "Displaying standard error ellipses") +
    theme_classic()+
    scale_color_discrete(labels = c(
      "high_high" = "High fibre + high exercise",
      "high_low"  = "High fibre + low exercise",
      "low_high"  = "Low fibre + high exercise",
      "low_low"   = "Low fibre + low exercise"
    ))
  
  ggsave(plot = p, filename = "results/aim1/06-ordination_plot.png", width = 8, height = 5)
  
  ############ LASSO
  # correlation heatmap
  meta_cor <- meta |>
    select(-fibre_group, -exercise_group, -fibre_exercise_group, -BMI_class, -Cardiometabolic_status, -city, -medication, -sex, -smoker, -stool_consistency)

  cor_matrix <- cor(meta_cor, use = "pairwise.complete.obs", method = "spearman")
  
  p_matrix <- matrix(NA, nrow = ncol(meta_cor), ncol = ncol(meta_cor))
  rownames(p_matrix) <- colnames(meta_cor)
  colnames(p_matrix) <- colnames(meta_cor)
  
  for (i in 1:ncol(meta_cor)) {
    for (j in 1:ncol(meta_cor)) {
      p_matrix[i, j] <- get_p(meta_cor[, i], meta_cor[, j])
    }
  }
  
  stars_matrix <- ifelse(p_matrix <= 0.001, "***",
                         ifelse(p_matrix <= 0.05, "**",
                                ifelse(p_matrix <= 0.1, "*",
                                       "")))
  

  png("results/aim1/0-heatmap_num_variables.png",
      width = 11.75, height = 9.52, units = "in", res = 300)
  
  pheatmap(
    cor_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = stars_matrix,
    number_color = "black",
    fontsize_number = 12,
    main = "Heatmap of Metadata Numeric Variables (Spearman correlation)"
  )
  
  dev.off()
  
  # mutating cv to "0" or "1". select numeric variables only
  meta_logistic <- meta |>
    mutate(Cardiometabolic_status = ifelse(Cardiometabolic_status == "Healthy", 0, 1))|>
    select(-fibre_group, -exercise_group, -fibre_exercise_group, -BMI_class, -city, -medication, -sex, -smoker, -stool_consistency)
  
  # splitting the data into a model selection and an inference set
  meta_split <- initial_split(meta_logistic, prop = 0.5, strata = Cardiometabolic_status)
  meta_training <- training(meta_split)
  meta_testing <- testing(meta_split)
  
  # model selection:
  
  # preparing the X and Y:
  X_training <-  model.matrix(Cardiometabolic_status ~  ., data=meta_training)[,-1]
  Y_training <- meta_training[,"Cardiometabolic_status"]=="1"
  
  # performing LASSO cross-validation to find lambda.min:
  meta_model <- 
    cv.glmnet(
      x = X_training, y = Y_training,
      alpha = 1,
      family = "binomial",
      type.measure = "auc",
      nfolds = 10)
  
  png("results/aim1/07-lasso_crossvalidation.png", width = 10, height = 10, units = "in", res = 300)
  plot(meta_model, main = "Cross-Validation with LASSO\n\n")
  dev.off()
  
  # meta_model$cvm[meta_model$index[1]]
  lambda.min <- meta_model$lambda.min
  
  # finding the covariates in lambda.min:
  lasso_model <- glmnet(
    x = X_training, y = Y_training,
    alpha = 1,
    family = "binomial",
    lambda = lambda.min
  )
  
  # obtaining the coefficients
  beta_lasso <- coef(lasso_model, s = "lambda.min")
  
  lasso_selected_covariates <- as_tibble(
    as.matrix(beta_lasso),
    rownames='covariate') |>
    filter(covariate != '(Intercept)' & abs(`s=lambda.min`) !=0)|>
    pull(covariate)
  
  # inference
  # fitting a logistic regression model with the variables selected by LASSO:
  lasso_formula <- reformulate(lasso_selected_covariates, response = "Cardiometabolic_status")
  
  lasso_log <- 
    glm(
      formula = lasso_formula,
      data = meta_testing,
      family = binomial)
  
  # hypothesis testing:
  lasso_log_results <- tidy(lasso_log, exponentiate = TRUE, conf.int = TRUE) |>
    mutate_if(is.numeric, round, 3)
  
  write_tsv(lasso_log_results, "results/aim1/08-lasso_hypothesis_testing.tsv")
  
  # prediction
  # create new model matrix for prediction:
  model_matrix_X_test <- 
    model.matrix(object = Cardiometabolic_status ~ .,
                 data = meta_testing)[, -1]
  
  # predict CV status on the test set using LASSO model:
  cv_predictions <- 
    round(predict(lasso_model, newx = model_matrix_X_test, type = "response"), 0)
  
  
  # obtain a confusion matrix for evaluation. "abnormal" as positive:
  confusion_matrix <- 
    confusionMatrix(
      data = as.factor(cv_predictions),
      reference = as.factor(meta_testing$Cardiometabolic_status),
      positive = "1"
    )
  
  saveRDS(confusion_matrix, "results/aim1/09-lasso_confusion_matrix.rds")
  
  # get ROC 
  ROC_lasso <- 
    roc(
      response = meta_testing$Cardiometabolic_status,
      predictor = predict(lasso_model,
                          newx = model_matrix_X_test)[,"s0"]) 
  
  png("results/aim1/10-lasso_roc_curve.png", width = 10, height = 10, units = "in", res = 300)
  plot(ROC_lasso, main = "ROC Curve for Selected LASSO Model", print.auc = TRUE)
  dev.off()
}

# helper functions
get_p <- function(x, y) {
  tmp <- cor.test(x, y, method = "spearman", exact = FALSE)
  tmp$p.value
}


main()

