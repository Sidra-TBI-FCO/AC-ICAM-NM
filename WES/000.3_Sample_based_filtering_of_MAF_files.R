
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("GenVisR", "dplyr", "ggplot2", "maftools", "ComplexHeatmap",
                      "data.table", "tidyr", "stringr", "cowplot", "tidyverse")
ipak(required.packages)

# Set parameters


# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_all_samples.Rdata")
load("./Processed_Data/WES/MAF/finalMaf_all_samples.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, pad = 0, 3)

# Keep only tumor samples
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Tissue == "T"),]
finalMaf = finalMaf[which(finalMaf$Tissue == "T"),]

# Exclude QC samples and non-epithelial
finalMafFiltered = finalMafFiltered[-which(finalMafFiltered$Patient_ID %in% excluded_df$Patient_ID),]
finalMaf = finalMaf[-which(finalMaf$Patient_ID %in% excluded_df$Patient_ID),]

save(finalMafFiltered, file = "./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")
save(finalMaf, file =  "./Processed_Data/WES/MAF/finalMaf_clean_primary_tumor_samples.Rdata")

# Exclude QC samples only for neoantigen prediction
load("./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_all_samples.Rdata")
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Tissue == "T"),]
finalMafFiltered = finalMafFiltered[-which(finalMafFiltered$Patient_ID %in% excluded_df$Patient_ID),]

save(finalMafFiltered, file = "./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_clean_samples.Rdata")
