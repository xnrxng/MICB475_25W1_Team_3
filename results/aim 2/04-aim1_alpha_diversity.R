# This script performs alpha diversity for aim 2.
# Date: November 15th 2025
# Usage: Rscript R/04-aim1_alpha_diversity.R

# Load the libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(picante)
library(ggplot2)
library(cowplot)
set.seed(2025)

# Read data
meta <- read.delim("data/data_raw/colombia_metadata.txt", sep = "\t")
otu <- read.delim(file = "data/data_processed/feature-table.txt", sep = "\t", skip = 1, row.names = 1)
taxa <- read.delim(file = "data/data_processed/taxonomy.tsv", sep = "\t")
tree <- read.tree(file = "data/data_processed/tree.nwk")

# Convert into df. Change rownames to be the sample names. Get rid of sample_name
meta_df <- as.data.frame(meta[,-1])
rownames(meta_df) <- meta$X.SampleID
sample <- sample_data(meta_df)

# Convert it into matrix 
otu_mat <- as.matrix(otu)
otu_tb <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Get rid of confidence column. Separate taxon into 7 columns. Convert into matrix
taxa_mat <- taxa |> 
  dplyr::select(-Confidence) |>
  separate(col = Taxon, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) |>
  as.matrix()

# Get rid of feature id as a column and make it the rownames
taxa_mat <- taxa_mat[,-1]
rownames(taxa_mat) <- taxa$`Feature.ID`
taxa_tb <- tax_table(taxa_mat)

# Convert it into the phyloseq object
phyloseq_obj <- phyloseq(otu_tb, sample, taxa_tb, tree)

# Keep bacteria. Get rid of chloroplasts and mitochondria
phyloseq_obj <- subset_taxa(phyloseq_obj, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Extract required components
OTU <- as.data.frame(otu_table(phyloseq_obj))
TREE <- phy_tree(phyloseq_obj)

# Calculate Faith’s PD
faith_pd <- picante::pd(OTU, TREE, include.root = TRUE)
