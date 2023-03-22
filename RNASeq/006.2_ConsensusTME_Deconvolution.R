
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

#install.packages("devtools")
#devtools::install_github("cansysbio/ConsensusTME")
library(ConsensusTME)

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load(paste0(toolbox.path, "/GSEA list/immune.gene.lists.v4.Rdata"))
source(paste0(toolbox.path, "/R scripts/consensusTMEAnalysis_WJed.R"))

# Testing package functions
#rawMethodSignatures = ConsensusTME::methodSignatures
#consensusGeneSets = ConsensusTME::consensusGeneSets

# Perform ssGSEA 
ES = consensusTMEAnalysis(RNASeq.QN.LOG2, cancer = "COAD", statMethod = "ssgsea")
save(ES, file = paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/ConsensusTME_COAD_ES.Rdata"))

#  Perform ssGSEA edited version
ES = consensusTMEAnalysis_WJed(RNASeq.QN.LOG2, cancer = "COAD", statMethod = "ssgsea", include_source = TRUE) # include_source = TRUE (when you want to see the source of the signature behind it)
save(ES, file = paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/ConsensusTME_COAD_WJedit_source_ES.Rdata"))
