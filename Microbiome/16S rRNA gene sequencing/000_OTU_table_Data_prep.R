

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

#if(!requireNamespace("BiocManager")){
# install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
#BiocManager::install("microbiome")
required.packages = c("stringr", "ggplot2", "phyloseq", "microbiome")
ipak(required.packages)
library(microbiome)

# load data
load("./Processed_Data/Microbiome/NormObjects.RData")
load("./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.2_sil_large_df.Rdata")

# Generate translation table
datastruct$SampleCode_new = paste(datastruct$SampleCode, datastruct$SamplesDescription, sep = "")
datastruct_small = datastruct[, c("Sample", "SampleCode_new")]
datastruct_small = datastruct_small[-which(duplicated(datastruct_small$Sample)),]
translation_table = datastruct_small
rownames(translation_table) = NULL
translation_table$SampleCode_new = str_pad(translation_table$SampleCode_new,
                                           4, pad = "0")
translation_table$SampleCode_new = gsub("LM1LM1", "LM1", translation_table$SampleCode_new)
translation_table$SampleCode_new = gsub("LM2LM2", "LM2", translation_table$SampleCode_new)
translation_table$SampleCode_new = gsub("291LM2LM1", "291LM2", translation_table$SampleCode_new)
translation_table = translation_table[order(translation_table$SampleCode_new),]

save(translation_table, file = "./Processed_Data/Microbiome/001_data_preparation/translation_table_Sample_Codes.Rdata")

# Summarizing the contents of a phyloseq object
summarize_phyloseq(s16sV1V3.2_sil)

# Compositional = NO
#1] Min. number of reads = 2869 
#2] Max. number of reads = 2869 
#3] Total number of reads = 1635330 
#4] Average number of reads = 2869 
#5] Median number of reads = 2869 
#7] Sparsity = 0.988239324847179 
#6] Any OTU sum to 1 or less? YES 
#8] Number of singletons = 6134 
#9] Percent of OTUs that are singletons 37.1442412498486 
#10] Number of sample variables are: 34 

# Generate OTU tables
OTU_table = otu_table(s16sV1V3.2_sil)
tax_table = tax_table(s16sV1V3.2_sil)
sample_data = sample_data(s16sV1V3.2_sil)
phy_tree(s16sV1V3.2_sil)
# Phylogenetic tree with 16514 tips and 16104 internal nodes.

# OTU (operational taxonomic units) table
# taxa are rows
OTU_table_dat = OTU_table@.Data
dim(OTU_table_dat)
# 16514   570

# Taxonomy Table: taxa by 7 taxonomic ranks
tax_table_dat = tax_table@.Data
dim(tax_table_dat)
# 16514     7
colnames(tax_table_dat)
# "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

sample_data_dat = data.frame(sample_data)

# Translate to proper sample ids
colnames(OTU_table_dat)

colnames(OTU_table_dat) = translation_table$SampleCode_new[match(colnames(OTU_table_dat), translation_table$Sample)]
OTU_table_dat = OTU_table_dat[, order(colnames(OTU_table_dat))]

dir.create("./Processed_Data/Microbiome/001_data_preparation/OTU_tables", showWarnings = FALSE)
dir.create("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Normalized_abundancies", showWarnings = FALSE)

save(OTU_table_dat, tax_table_dat, sample_data_dat, file = "./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Normalized_abundancies/s16sV1V3.2_sil_OTU_table.Rdata")

myTaxa = names(sort(taxa_sums(s16sV1V3.2_sil), decreasing = TRUE)[1:50])
ex1 = prune_taxa(myTaxa, s16sV1V3.2_sil)
plot(phy_tree(ex1), show.node.label = TRUE, show.tip.label = FALSE, color = "SamplesDescription",
     label.tips = "Phylum", ladderize = "left", justify = "left" , size = "Abundance")
