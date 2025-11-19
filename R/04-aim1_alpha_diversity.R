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
library(phangorn)
library(phytools)
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

# Resolve polytomies (make the tree dichotomous)
TREE_dich <- ape::multi2di(TREE)
# Check if branch lengths exist
head(TREE_dich$edge.length)
# Force ultrametric
ape::is.ultrametric(TREE_dich)
TREE_ultra <- ape::compute.brlen(TREE_dich, method = "Grafen")

# Ensure tree is rooted
TREE_rooted <- ape::root(TREE_dich, outgroup = TREE_dich$tip.label[1], resolve.root = TRUE)
# Make it fully bifurcating
TREE_rooted_dich <- ape::multi2di(TREE_rooted)
# Use phytools::chronos() clone = chronopl (so install + use phytools)
TREE_ultra <- phytools::chronopl(TREE_rooted_dich, lambda = 1)

# Ensure tree is dichotomous
TREE_dich <- ape::multi2di(TREE_rooted_dich)
# Apply Grafen branch lengths
TREE_ultra <- ape::compute.brlen(TREE_dich, method = "Grafen")
# Check
ape::is.ultrametric(TREE_ultra)
# Fix zero-length branches
TREE_ultra$edge.length[TREE_ultra$edge.length == 0] <- 1e-6

# Compute node ages using branching.times()
node_ages <- ape::branching.times(TREE_ultra)
# Attach ages to the tree for picante
ntip  <- length(TREE_ultra$tip.label)
nnode <- TREE_ultra$Nnode
TREE_ultra$ages <- c(
  rep(0, ntip),    # Tips have age 0
  node_ages        # Internal nodes get branching times
)
names(TREE_ultra$ages) <- 1:(ntip + nnode)

# Fix the OTU table orientation
OTU_pd <- as.data.frame(t(OTU))
# Check that columns = taxa
all(colnames(OTU_pd) == TREE_ultra$tip.label)
# Calculate Faith’s PD
faith_pd <- picante::pd(OTU_pd, TREE_ultra, include.root = TRUE)

