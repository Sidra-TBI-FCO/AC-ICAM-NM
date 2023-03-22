
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Load data
merged_txt = read.csv("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/8_merge_maf_files/Mutect/samples.maf", sep = "\t")
head.txt = read.table("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/8_merge_maf_files/Mutect/head_samples.txt", sep = "\t", stringsAsFactors = FALSE)
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv")

# Set parameters
version = "3"

# Set colnames
colnames(merged_txt) = head.txt[1,]

MAF_df = merged_txt

# Add Patient ID column
MAF_df$Sample_ID = gsub("WES_COAD_LUMC_SIDRA_", "", MAF_df$Tumor_Sample_Barcode)
MAF_df$Sample_ID = gsub(".sorted", "", MAF_df$Sample_ID)
MAF_df$Patient_ID = gsub("T", "", MAF_df$Sample_ID)
MAF_df$Patient_ID = gsub("LM1", "", MAF_df$Patient_ID)
MAF_df$Patient_ID = gsub("LM2", "", MAF_df$Patient_ID)
MAF_df$Patient_ID = gsub("LM", "", MAF_df$Patient_ID)
MAF_df$Tissue = str_sub(MAF_df$Sample_ID,-1,-1)
MAF_df$Tissue[which(MAF_df$Tissue == "M")] = "LM"
MAF_df$Tissue[which(MAF_df$Tissue == "1")] = "LM1"
MAF_df$Tissue[which(MAF_df$Tissue == "2")] = "LM2"

MAF_df$Patient_ID = str_pad(MAF_df$Patient_ID, 3, pad = "0") # Add leading zeros

MAF_df$Sample_ID = paste(MAF_df$Patient_ID, MAF_df$Tissue, sep = "")
MAF_df$FA = NULL

mutect_MAF_df = MAF_df

# Save as R data file
dir.create("./Processed_Data/WES/MAF", showWarnings = FALSE)
#save(MAF_df, file = paste0("./Processed_Data/WES/MAF/MAF_Primary_tumors_unfiltered_version", version, ".Rdata"))
save(mutect_MAF_df, file = paste0("./Processed_Data/WES/MAF/MAF_mutect_all_version_", version, ".Rdata"))

