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
}

#helper functions
main()