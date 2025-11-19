library(tidyverse)
library(tidymodels)
library(phyloseq)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)

set.seed(2025)

#read data
phyloseq <- readRDS("data/data_raw/phyloseq_object.rds")

##### Phyloseq Preprocessing #####

# change the rank names from "1,2.." to "Kingdom, Phlyum..." instead of "Rank 1, Rank 2..."
colnames(tax_table(phyloseq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(phyloseq) 


#Convert fibre intake and METs per week into categorical variables
sample_data(phyloseq)$Fibre_cat <- ifelse(sample_data(phyloseq)$fiber > 20,
                                        "Adequate", "Inadequate")

sample_data(phyloseq)$Exercise_cat <- ifelse(sample_data(phyloseq)$MET_mins_per_week > 1000,
                                           "Adequate", "Inadequate")

#Convert to factor
sample_data(phyloseq)$Fibre_cat    <- factor(sample_data(phyloseq)$Fibre_cat,
                                           levels = c("Inadequate","Adequate"))

sample_data(phyloseq)$Exercise_cat <- factor(sample_data(phyloseq)$Exercise_cat,
                                           levels = c("Inadequate","Adequate"))

#Confirm conversion
table(sample_data(phyloseq)$Fibre_cat)
table(sample_data(phyloseq)$Exercise_cat)

#Remove any ASV seen less at <0.1% abundance across the entire study
x = taxa_sums(phyloseq)
keepTaxa = (x / sum(x)) > 0.001 # 0.001 = 0.1%. Change this if you want a different cutoff. 
phyloseq_pruned = prune_taxa(keepTaxa, phyloseq)
#Check the resulting new phyloseq object
phyloseq_pruned

#Subset by CV status
phyloseq_CV <- subset_samples(phyloseq_pruned, Cardiometabolic_status %in% c("Abnormal","Healthy"))
phyloseq_CV <- prune_taxa(taxa_sums(phyloseq_CV) > 0, phyloseq_CV)

#Rename column in the subsetted phyloseq object
sample_data(phyloseq_CV)$CV_Status <- sample_data(phyloseq_CV)$Cardiometabolic_status

#Convert to factor
sample_data(phyloseq_CV)$CV_Status <- factor(sample_data(phyloseq_CV)$CV_Status,
                                             levels = c("Healthy","Abnormal"))

#Add 1 to all counts because DESeq does not work well with zeros
otu_table(phyloseq_CV) <- otu_table(otu_table(phyloseq_CV) + 1, taxa_are_rows = FALSE)

##### DESeq #####

#Convert to DESeq dataset
#Microbiome vs CV status, taking into account fibre and exercise
deseq_cv <- phyloseq_to_deseq2(
  phyloseq_CV,
  ~ Fibre_cat + Exercise_cat + CV_Status
)

#Run the DESeq2 model
deseq_cv <- DESeq(deseq_cv)

#Extract results for CV_Status Abnormal vs Healthy
res_cv <- results(deseq_cv, contrast = c("CV_Status","Abnormal","Healthy"))
res <- as.data.frame(res_cv)

#Extract genus labels instead of just ASV IDs
taxdf <- as.data.frame(tax_table(phyloseq_CV))

# ensure rownames in both objects match ASV IDs
asv_ids <- rownames(res)
taxdf_sub <- taxdf[asv_ids, , drop = FALSE]

#Create human readable labels
make_label <- function(x){
  # use Genus + Species if available
  if(!is.na(x["Genus"]) && x["Genus"] != ""){
    if(!is.na(x["Species"]) && x["Species"] != ""){
      return(paste(x["Genus"], x["Species"]))
    } else {
      return(x["Genus"])
    }
  }
  # otherwise fall back to the deepest available taxonomic rank
  deepest <- tail(na.omit(x), 1)
  if(length(deepest) > 0) return(as.character(deepest))
  
  return(NA_character_) 
}

res$Label <- apply(taxdf_sub, 1, make_label)

# fallback to ASV ID if still NA
res$Label[is.na(res$Label)] <- rownames(res)[is.na(res$Label)
                                             

##### Plots #####                                             
                                             
#Generating the volcano plot
EnhancedVolcano(as.data.frame(res),
                lab = res$Label, 
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Abnormal vs Healthy (CV Status)")


##Generating bar plot
ggplot(res, aes(x = reorder(Label, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_manual(values = c("red", "blue"), labels = c("Downregulated","Upregulated")) +
  labs(title = "Differentially Abundant Taxa in Abnormal CV Status",
       x = "Taxa",
       y = "log2 Fold Change (Abnormal vs Healthy)",
       fill = "Direction") +
  theme_minimal(base_size = 8)
