
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/stefanofunctions.R"))                                                                 # Used for calinsky function and plot

# Set Parameters

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")

## Subset to remove liver meta samples and biological replicates
liver_meta_samples = colnames(RNASeq.QN.LOG2)[grep("LM", colnames(RNASeq.QN.LOG2))]
biological_replicates = colnames(RNASeq.QN.LOG2)[c(grep("B", colnames(RNASeq.QN.LOG2)), grep("C", colnames(RNASeq.QN.LOG2)))]
to_remove = c(liver_meta_samples, biological_replicates)
RNASeq.QN.LOG2 = RNASeq.QN.LOG2[, -which(colnames(RNASeq.QN.LOG2) %in% to_remove)]
RNASeq.NORM.quantiles = RNASeq.NORM.quantiles[,-which(colnames(RNASeq.NORM.quantiles) %in% to_remove)]

save(RNASeq.QN.LOG2, RNASeq.NORM.quantiles, file = "./Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")

