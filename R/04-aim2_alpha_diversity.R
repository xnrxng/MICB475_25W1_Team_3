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
library(dplyr)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(ggsignif)
set.seed(2025)

main <- function(){}
# Read data
phyloseq_obj <- readRDS("data/data_processed/phyloseq_obj_rarefied.rds")

# Extract required components
OTU <- as.data.frame(otu_table(phyloseq_obj))
TREE <- phy_tree(phyloseq_obj)

# Resolve polytomies
TREE_dich <- ape::multi2di(TREE)

# Apply Grafen branch lengths (ultrametric)
TREE_ultra <- ape::compute.brlen(TREE_dich, method = "Grafen")

# Fix any zero-length edges
TREE_ultra$edge.length[TREE_ultra$edge.length == 0] <- 1e-6

# Calculate node ages if running picante functions requiring ages
node_ages <- ape::branching.times(TREE_ultra)
ntip  <- length(TREE_ultra$tip.label)
nnode <- TREE_ultra$Nnode

TREE_ultra$ages <- c(
  rep(0, ntip),
  node_ages
)
names(TREE_ultra$ages) <- 1:(ntip + nnode)

# Fix the OTU table orientation
OTU_pd <- as.data.frame(t(OTU))
# Check that columns = taxa
# all(colnames(OTU_pd) == TREE_ultra$tip.label)
# Calculate Faith’s PD (data frame (faith_pd) with Faith’s PD for each sample (rows correspond to samples, columns include PD and SR: species richness))
faith_pd <- picante::pd(OTU_pd, TREE_ultra, include.root = TRUE)

# Merge faith_pd with metadata
meta <- sample_data(phyloseq_obj)
meta$SampleID <- rownames(meta)
meta <- data.frame(meta)
faith_pd$SampleID <- rownames(faith_pd)
meta_pd <- merge(meta, faith_pd, by = "SampleID")
# nrow(meta_pd)  # Should now match number of samples

# Create the lifestyle_group column
meta_pd <- meta_pd |>
  mutate(
    lifestyle_group = case_when(
      fiber >= 20 & MET_mins_per_week > 1000 ~ "adequate fibre;high exercise",
      fiber >= 20 & MET_mins_per_week <= 1000  ~ "adequate fibre;low exercise",
      fiber < 20  & MET_mins_per_week > 1000 ~ "inadequate fibre;high exercise",
      TRUE ~ "inadequate fibre;low exercise"
    )
  )
# Check if it worked
# table(meta_pd$lifestyle_group)
# Make sure CV_status exists
# table(meta_pd$Cardiometabolic_status)
# Create the 8-group factor
meta_pd$group <- paste(meta_pd$lifestyle_group, meta_pd$Cardiometabolic_status, sep = "_")
# table(meta_pd$group)  # Check counts for each of 8 groups

# Save data
write_tsv(meta_pd, "results/aim2/alpha_diversity/01-meta_pd_faith_pd.tsv")
write_tsv(faith_pd, "results/aim2/alpha_diversity/02-faith_pd.tsv")

# Kruskal-Wallis test
healthy_pd <- meta_pd |>
  filter(Cardiometabolic_status == "Healthy")
kw_res_healthy <- kruskal.test(PD ~ group, data = healthy_pd)

abnormal_pd <- meta_pd |>
  filter(Cardiometabolic_status == "Abnormal")
kw_res_abnormal <- kruskal.test(PD ~ group, data = abnormal_pd)

# Save a tidy summary
kw_summary_healthy <- data.frame(
  statistic = kw_res_healthy$statistic,
  df = kw_res_healthy$parameter,
  p_value = kw_res_healthy$p.value
)

kw_summary_abnormal <- data.frame(
  statistic = kw_res_abnormal$statistic,
  df = kw_res_abnormal$parameter,
  p_value = kw_res_abnormal$p.value
)

kw_summary <- rbind(kw_summary_healthy, kw_summary_abnormal)
rownames(kw_summary) <- c("healthy", "abnormal")

write_tsv(kw_summary, "results/aim2/alpha_diversity/03-kruskal_wallis_group.tsv")

# More saving
pairwise_res <- pairwise.wilcox.test(meta_pd$PD, meta_pd$group, p.adjust.method = "BH")
saveRDS(pairwise_res, "results/aim2/alpha_diversity/04-pairwise_wilcox_group.rds")
# Tidy version for easier reading
pairwise_df <- as.data.frame(pairwise_res$p.value) |>
  tibble::rownames_to_column("Group1")
write_tsv(pairwise_df, "results/aim2/alpha_diversity/05-pairwise_wilcox_group.tsv")

# Boxplot
kw_summary <- kw_summary |>
  mutate(label = paste0("KW p = ", signif(p_value, 3)))

kw_summary$Cardiometabolic_status <- c("Healthy", "Abnormal")

bp <- ggplot(meta_pd, aes(x = lifestyle_group, y = PD, fill = lifestyle_group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  facet_wrap(~ Cardiometabolic_status, scales = "free_x") +
  labs(
    x = NULL,
    y = "Faith's PD",
    title = "Faith's PD by Lifestyle and CV Status"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)) +
  scale_fill_brewer(
    palette = "Set2")+
  scale_x_discrete(labels = c(
    "adequate fibre;high exercise"    = "High fibre\nHigh exercise",
    "adequate fibre;low exercise"     = "High fibre\nLow exercise",
    "inadequate fibre;high exercise"  = "Low fibre\nHigh exercise",
    "inadequate fibre;low exercise"   = "Low fibre\nLow exercise"
  )) +
  geom_text(
    data = kw_summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) 
# Save
ggsave("results/aim2/alpha_diversity/06-faith_PD_boxplot.png", 
       plot = bp, 
       width = 10, height = 6, units = "in", dpi = 300)
}

# Calculate Shannon's
shannon_df <- phyloseq::estimate_richness(phyloseq_obj, measures = "Shannon")
# Estimate_richness returns a data.frame with rownames = sample IDs, column "Shannon"
shannon_df <- tibble::rownames_to_column(shannon_df, var = "SampleID")
# Keep only the Shannon column
shannon_df <- shannon_df[, c("SampleID", "Shannon")]

# Ensure SampleID columns are character
meta_pd$SampleID <- as.character(meta_pd$SampleID)
shannon_df$SampleID <- as.character(shannon_df$SampleID)

# Merge (inner join)
meta_pd <- dplyr::left_join(meta_pd, shannon_df, by = "SampleID")

# Check
dplyr::glimpse(meta_pd)
table(is.na(meta_pd$Shannon))   # should be 0 FALSE ideally

# Kruskal-Wallis test
healthy_pd <- meta_pd |>
  filter(Cardiometabolic_status == "Healthy")
kw_shannon_healthy <- kruskal.test(Shannon ~ group, data = healthy_pd)

abnormal_pd <- meta_pd |>
  filter(Cardiometabolic_status == "Abnormal")
kw_shannon_abnormal <- kruskal.test(Shannon ~ group, data = abnormal_pd)

# Save a tidy summary
kw_summary_healthy2 <- data.frame(
  statistic = kw_shannon_healthy$statistic,
  df = kw_shannon_healthy$parameter,
  p_value = kw_shannon_healthy$p.value
)

kw_summary_abnormal2 <- data.frame(
  statistic = kw_shannon_abnormal$statistic,
  df = kw_shannon_abnormal$parameter,
  p_value = kw_shannon_abnormal$p.value
)

kw_summary2 <- rbind(kw_summary_healthy2, kw_summary_abnormal2)
rownames(kw_summary2) <- c("healthy", "abnormal")
write_tsv(kw_summary2, "results/aim2/alpha_diversity/07-kruskal_wallis_group2.tsv")

# More saving
pairwise_shannon <- pairwise.wilcox.test(meta_pd$Shannon, meta_pd$group, p.adjust.method = "BH")
saveRDS(pairwise_shannon, "results/aim2/alpha_diversity/08-pairwise_wilcox_group2.rds")
# Tidy version for easier reading
pairwise_shannon_df <- as.data.frame(pairwise_shannon$p.value) |>
  tibble::rownames_to_column("Group1")
write_tsv(pairwise_shannon_df, "results/aim2/alpha_diversity/09-pairwise_wilcox_group2.tsv")

# Boxplot
bp2 <- ggplot(meta_pd, aes(x = lifestyle_group, y = Shannon, fill = lifestyle_group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  facet_wrap(~ Cardiometabolic_status, scales = "free_x") +
  labs(
    x = NULL,
    y = "Shannon diversity (H')",
    title = "Shannon diversity by Lifestyle and CV Status"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.spacing = unit(0, "lines"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)) +
  scale_fill_brewer(
    palette = "Set2")+
  scale_x_discrete(labels = c(
    "adequate fibre;high exercise"    = "High fibre\nHigh exercise",
    "adequate fibre;low exercise"     = "High fibre\nLow exercise",
    "inadequate fibre;high exercise"  = "Low fibre\nHigh exercise",
    "inadequate fibre;low exercise"   = "Low fibre\nLow exercise"
  )) +
  geom_text(
    data = kw_summary2,
    aes(x = -Inf, y = Inf, label = p_value),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) 

# Save
ggsave("results/aim2/alpha_diversity/10-Shannon_boxplot.png", 
       plot = bp, 
       width = 10, height = 6, units = "in", dpi = 300)
}

main()
