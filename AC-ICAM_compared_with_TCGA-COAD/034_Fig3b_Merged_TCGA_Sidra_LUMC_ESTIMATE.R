
# Run ESTIMATE algorithm on merged dataset (TCGA + Sidra-LUMC)

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)

library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(estimate)

# Load data
load("./Processed_Data/RNASeq/030_Merged_TCGA_COAD_Sidra_LUMC/TCGA_COAD_Sidra_LUMC_Primary_tumor_EDASEQ_QN_LOG2.Rdata")
write.table(RNASeq.QN.LOG2, sep = "\t", file = "./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_merged_expression_data.txt", quote = FALSE)

dir.create("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU", showWarnings = FALSE)


# Calculate estimate score for Merged cohort
filterCommonGenes(input.f= paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_merged_expression_data.txt"),
                  output.f= paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.input.gct"),
                  id=c("GeneSymbol","EntrezID"))

estimateScore(input.ds = paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.input.gct"),
              output.ds = paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.score.gct"),
              platform= "illumina")

estimate.gct<-read.table(paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.score.gct"), skip = 2, header = TRUE) #skip=2 remove the first 2 lines
rownames(estimate.gct) = estimate.gct$NAME
estimate.gct$NAME = NULL
estimate.gct$Description = NULL

ESTIMATE = t(estimate.gct)
rownames(ESTIMATE) = gsub("X", "", rownames(ESTIMATE))
rownames(ESTIMATE) = gsub("\\.", "-", rownames(ESTIMATE))
save(ESTIMATE, file = paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.score.Rdata"))
