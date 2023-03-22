
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

#if(!requireNamespace("BiocManager")){
 # install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
required.packages = c("stringr", "ggplot2", "phyloseq")
ipak(required.packages)

# load data
load("./Processed_Data/Microbiome/NormObjects.RData")

# Use phyloseq psmelt function to convert to datastructures
datastruct = psmelt(s16sV1V3.2_sil)

dir.create("./Analysis/Microbiome", showWarnings = FALSE)
dir.create("./Processed_Data/Microbiome/001_data_preparation", showWarnings = FALSE)

save(datastruct, file = "./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.2_sil_large_df.Rdata")

# Use phyloseq psmelt function to convert to datastructures
datastruct_ICR = psmelt(s16sV1V3.4_sil_icr)
save(datastruct_ICR, file = "./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.4_sil_icr_large_df.Rdata")

datastruct_ajcc = psmelt(s16sV1V3.4_sil_ajcc)
save(datastruct_ajcc, file = "./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.4_sil_ajcc_large_df.Rdata")

datastruct_tissue = psmelt(s16sV1V3.4_sil_tissue)
save(datastruct_tissue, file = "./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.4_sil_tissue_large_df.Rdata")
