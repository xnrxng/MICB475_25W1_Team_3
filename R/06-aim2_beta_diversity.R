# This script performs beta diversity under aim 2. 
# Date: November 16, 2025
# Usage: Rscript R/06-aim2_beta_diversity.R

# Load in libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)
library(picante)
library(ggplot2)
library(cowplot)
library(viridis)
library(tidyr) #For some reason the 'separate' function wouldn't work unless I loaded tidyr in separately (haha, nice)

# Set the seed for reproducibility 
set.seed(2025)

# Read in data
meta <- read.delim("data/data_raw/colombia_metadata.txt", sep = "\t")
otu <- read.delim(file = "data/data_processed/feature-table.txt", sep = "\t", skip = 1, row.names = 1)
taxa <- read.delim(file = "data/data_processed/taxonomy.tsv", sep = "\t")
tree <- read.tree(file = "data/data_processed/tree.nwk")

# Prepare metadata by converting to a df. Use SampleID as row name. Remove SampleID from first column.
meta_df <- meta %>% select(-X.SampleID)
rownames(meta_df) <- meta$X.SampleID

# Create the lifestyle_group (The 4 categories)
# NOTE the metadata spells "Fibre" as "Fiber" so be careful about spelling when calling on the category.
meta_lg <- meta_df %>%
  mutate(lifestyle_group = case_when(
    fiber >= 20 & MET_mins_per_week >= 1000 ~ "adequate fibre;high exercise",
    fiber >= 20 & MET_mins_per_week < 1000 ~ "adequate fibre;low excercies",
    fiber < 20 & MET_mins_per_week >= 1000 ~ "inadequate fibre;high excercise",
    TRUE ~ "inadequate fibre;low exercise"
  ))

# Check if it worked
table(meta_lg$lifestyle_group)
# Double check to ensure that Cardiometabolic_status exists
table(meta_df$Cardiometabolic_status)

# Create the 8 group factor by subdivide each of the four composite groups by CV health (healthy;poor).  
meta_gf <- meta_lg %>%
  mutate(group = paste(lifestyle_group, Cardiometabolic_status, sep = "_"))

sample <- sample_data(meta_gf)

# Check counts for the 8 groups
table(meta_gf$group)

# Save metadata for later! Woohoo~ 
saveRDS(meta_gf, "results/aim2/beta_diversity/00-meta_grouped.rds")
write_tsv(meta_gf, "results/aim2/beta_diversity/01-meta_grouped.tsv")

# Prepare the OTU table
otu_mat <- as.matrix(otu)
otu_tb <- otu_table(otu_mat, taxa_are_rows = TRUE)


# Prepare Taxonomy Table: Remove confidence, separate taxon into 7 columns, convert to matrix. YIPEEEEE. P.S I know that these comments will get deleted but they make me happy so thanks for reading em'.
taxa_mat <- taxa %>%
  dplyr::select(-Confidence) %>%
  separate(col = Taxon, sep = "; ",
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()

# Get rid of feature id as a column and make it the rownames
taxa_mat <- taxa_mat[,-1]
rownames(taxa_mat) <- taxa$`Feature.ID`
taxa_tb <- tax_table(taxa_mat)
  

# Generate the phyloseq object
phylo_obj <- phyloseq(otu_tb, sample, taxa_tb, tree)

# Remove Chloroplasts/Mitochondria. Keep aaaallllllll the Bacteria.
phylo_obj <- subset_taxa(phylo_obj,Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# ------------------Gets Funky A Bit----------------------------------------

#PCoA ordination
bray_ord <- ordinate(phylo_obj, method = "PCoA", distance = "bray")

#plot
bray_plot <- plot_ordination(phylo_obj, bray_ord, color = "Cardiometabolic_status") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Bray-Curtis PCoA by CV Status",
       color = "Cardiometabolic_status")

print(bray_plot)

#----------------------------------------------------------------------------
### TRYING SUM - OKAY I THINK I COOKED ON THIS ONE DOWN BELOW***

#Prep OTU 
OTU <- as.data.frame(t(otu_table(phylo_obj)))

#Compute Bray-Curtis
bray_dist <- vegdist(OTU, method = "bray")

#PCoA ordination (Bray–Curtis)
bray_ord <- ordinate(phylo_obj, method = "PCoA", distance = "bray")

#Plot PCoA by 8 groups
plot_pcoa <- plot_ordination(
  phylo_obj,
  bray_ord,
  color = "group"
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_viridis_d(option = "turbo") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Bray–Curtis PCoA of Gut Microbiome",
    color = "Lifestyle × CV Status"
  )

print(plot_pcoa)

#Save the plot
ggsave("results/aim2/beta_diversity/02-Bray_Curtis_PCoA_Plot.png", 
       plot = plot_pcoa, 
       width = 10, height = 6, units = "in", dpi = 300)


#Test beta dispersion
# Extract grouping variable
group8 <- sample_data(phylo_obj)$group

# Beta dispersion test
beta_disp <- betadisper(bray_dist, group8)
anova(beta_disp)
# Save a lil copy
saveRDS(beta_disp, "results/aim2/beta_diversity/03-beta_dispersion.rds")

# PERMANOVA: Compare Bray–Curtis across 8 groups
permanova_results <- adonis2(
  bray_dist ~ group8,
  data = as(sample_data(phylo_obj), "data.frame"),
  permutations = 9999
)

print(permanova_results)
# Save a lil copy
saveRDS(permanova_results, "results/aim2/beta_diversity/04-permanova_results.rds")


