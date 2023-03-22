
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

#install.packages("devtools")
#devtools::install_github("cansysbio/ConsensusTME")
library(ConsensusTME)

# Load data
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)

# Testing package functions
#rawMethodSignatures = ConsensusTME::methodSignatures
#consensusGeneSets = ConsensusTME::consensusGeneSets

# Perform ssGSEA 
ES = consensusTMEAnalysis(RNASeq.QN.LOG2, cancer = "COAD", statMethod = "ssgsea")
save(ES, file = paste0("./Analysis/Deconvolution_and_GSEA/ConsensusTME_COAD_ES.Rdata"))
