# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "reshape2", "ComplexHeatmap", "circlize", "dendsort", "openxlsx", "dplyr"))

load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
model = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv")
load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_patients_microbiome.Rdata")

##########
# fix TCGA 
microbiome.TCGA = tcga_bacteria_df
rownames(microbiome.TCGA) = microbiome.TCGA$name
microbiome.TCGA$name = NULL

# colnames(microbiome.TCGA) = gsub("\\.", "-", colnames(microbiome.TCGA))
# colnames(microbiome.TCGA) = substring(colnames(microbiome.TCGA), 1, 12)

# edit modle genera names
model$covariates = gsub(".*\\D_5__", "", model$covariates)

#########
# add genera 
matrix.16s = as.data.frame(Genus_full_abundance)
matrix.16s$genera = gsub(".*\\D_5__", "", rownames(matrix.16s))

matrix.WGS = Genus_WGS
rownames(matrix.WGS) = gsub("g__", "", rownames(matrix.WGS))
rownames(matrix.WGS) = gsub("_", " ", rownames(matrix.WGS))

#########
# subset 16S and WGS
matrix.16s = matrix.16s[which(matrix.16s$genera %in% c(rownames(matrix.WGS), "Ruminococcus 1", "Ruminococcus 2")),]
matrix.WGS = matrix.WGS[which(rownames(matrix.WGS) %in% c(matrix.16s$genera, "Ruminococcus")),]

########
# subset TCGA 
microbiome.TCGA = microbiome.TCGA[-which(rowSums(microbiome.TCGA) == 0),]
microbiome.TCGA = microbiome.TCGA[which(rownames(microbiome.TCGA) %in% c(matrix.16s$genera, "Ruminococcus")),]

################################################################################
################################################################################

# overlap model genera 
model.16s = matrix.16s[which(matrix.16s$genera %in% model$covariates),]
model.wgs = matrix.WGS[which(rownames(matrix.WGS) %in% c(model$covariates, "Ruminococcus")),]


################################################################################
################################################################################

# TCGA overlap with 16S

rownames(tcga_bacteria_df) = tcga_bacteria_df$name
tcga_bacteria_df$name = NULL
tcga = tcga_bacteria_df[-which(rowSums(tcga_bacteria_df) == 0),]

genera.16s = as.data.frame(Genus_full_abundance)
genera.16s$genera = gsub(".*\\D_5__", "", rownames(genera.16s))
tcga = tcga[which(rownames(tcga) %in% c(genera.16s$genera, "Ruminococcus")),]

mode.tcga = tcga[which(rownames(tcga) %in% c(model$covariates, "Ruminococcus")),]
