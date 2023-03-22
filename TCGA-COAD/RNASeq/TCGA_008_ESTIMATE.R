#################################################################
###
### 
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(estimate)

# Set parameters
download.method = "TCGA_Assembler"

# Load data
if(download.method == "TCGA_Assembler"){
  load("./Processed_Data/RNASeq/TCGA_Assembler_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
}
if(download.method == "Biolinks"){
  load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
}

write.table(filtered.norm.RNAseqData, sep = "\t", file = paste0("./Processed_Data/RNASeq/", download.method, "_RNASeq_for_ESTIMATE", ".txt"), quote = FALSE)
dir.create("./Analysis/008_ESTIMATE", showWarnings = FALSE)

# Calculate estimate score for TCGA COAD cohort
filterCommonGenes(input.f=paste0("./Processed_Data/RNASeq/", download.method, "_RNASeq_for_ESTIMATE", ".txt"),
                  output.f=paste0("./Analysis/008_ESTIMATE/", download.method, ".input.gct"),
                  id=c("GeneSymbol","EntrezID"))

estimateScore(input.ds = paste0("./Analysis/008_ESTIMATE/", download.method, ".input.gct"),
              output.ds = paste0("./Analysis/008_ESTIMATE/", download.method, "_ESTIMATE_score.gct"),
              platform= "illumina")

estimate.gct<-read.table(paste0("./Analysis/008_ESTIMATE/", download.method, "_ESTIMATE_score.gct"), skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate.gct) = estimate.gct$NAME
estimate.gct$NAME = NULL
estimate.gct$Description = NULL

ESTIMATE = t(estimate.gct)
rownames(ESTIMATE) = gsub("\\.", "-", rownames(ESTIMATE))
save(ESTIMATE, file = paste0("./Analysis/008_ESTIMATE/", download.method, "_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata"))
