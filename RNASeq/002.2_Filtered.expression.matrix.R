
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/stefanofunctions.R"))  
ipak("stringr")

# Load
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Filter
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -which(substring(colnames(RNASeq.QN.LOG2), 1, 3) %in% excluded_df$Patient_ID)]
dim(RNASeq.QN.LOG2)
RNASeq.NORM.quantiles = RNASeq.NORM.quantiles[, -which(substring(colnames(RNASeq.NORM.quantiles), 1, 3) %in% excluded_df$Patient_ID)]
dim(RNASeq.NORM.quantiles)

# Save
save(RNASeq.NORM.quantiles, RNASeq.QN.LOG2, file = "./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")