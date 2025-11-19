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
library(microbiome)
set.seed(2025)

main <- function(){
  #read data
  phylo_rare <- readRDS("data/data_processed/phyloseq_obj_rarefied.rds")
  
  # aggregate
  ps_phylum <- tax_glom(phylo_rare, taxrank = "Phylum")
  
  # TSS transform
  ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
  
  # CLR transform
  ps_phylum_clr <- microbiome::transform(ps_phylum_rel, "clr")
  
  # prepare for lasso
  phylum_tb <- otu_table(ps_phylum_clr)
  tax_table_data <- tax_table(ps_phylum)
  rownames(phylum_tb) <- tax_table_data[, "Phylum"] |> as.vector()
  phylum_tb <- as.data.frame(phylum_tb)
  
  meta <- sample_data(ps_phylum)
  meta <- data.frame(meta)
  
  # make sure rows are in the same order
  phylum_tb <- t(phylum_tb)
  phylum_tb <- phylum_tb[rownames(meta), , drop = FALSE]
  
  # mutating cv to "0" or "1". append to phylum_tb
  meta_logistic <- meta |>
    mutate(Cardiometabolic_status = ifelse(Cardiometabolic_status == "Healthy", 0, 1)) |>
    select(Cardiometabolic_status)
  
  phylum_meta <- cbind(phylum_tb, meta_logistic)
  
  # splitting the data into a model selection and an inference set
  meta_split <- initial_split(phylum_meta, prop = 0.5, strata = Cardiometabolic_status)
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
  
  png("results/aim3/01-lasso_crossvalidation.png", width = 10, height = 10, units = "in", res = 300)
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
  
  write_tsv(lasso_log_results, "results/aim3/02-lasso_hypothesis_testing.tsv")
  
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
  
  saveRDS(confusion_matrix, "results/aim3/03-lasso_confusion_matrix.rds")
  
  # get ROC 
  ROC_lasso <- 
    roc(
      response = meta_testing$Cardiometabolic_status,
      predictor = predict(lasso_model,
                          newx = model_matrix_X_test)[,"s0"]) 
  
  png("results/aim3/4-lasso_roc_curve.png", width = 10, height = 10, units = "in", res = 300)
  plot(ROC_lasso, main = "ROC Curve for Selected LASSO Model", print.auc = TRUE)
  dev.off()
  
  
  
}

main()