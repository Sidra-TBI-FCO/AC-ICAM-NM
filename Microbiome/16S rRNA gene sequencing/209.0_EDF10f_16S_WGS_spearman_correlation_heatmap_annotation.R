# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# load packages
ipak(c("ggplot2", "ggpubr", "stringr", "corrplot", "openxlsx", "circlize", "ComplexHeatmap"))

# load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
tumor.RA = Genus_full_abundance
#WGS.filtered = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/WGS_validation/Genera_filtered_Relative_abundance_matrix_WGS_167.csv")
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
Best_CV_coef = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv", stringsAsFactors = FALSE)

# edit names 
Best_CV_coef$covariates = gsub(".*\\D_5__", "", Best_CV_coef$covariates)

# subset 16s to model genera
tumor.RA = tumor.RA[which(rownames(tumor.RA) %in% Best_CV_coef$coefficients),]
tumor.RA = as.data.frame(tumor.RA)
tumor.RA$genus = gsub(".*\\D_5__", "", rownames(tumor.RA))

# edit wgs genus names 
rownames(Genus_WGS) = gsub("g__", "", rownames(Genus_WGS))

# species 
Species_WGS = as.data.frame(t(Species_WGS))

# add bromii to match rum2
Genus_WGS = as.data.frame(t(Genus_WGS))
Genus_WGS$`Ruminococcus 2` = Species_WGS$s__Ruminococcus_bromii[match(rownames(Genus_WGS), rownames(Species_WGS))]
Genus_WGS$`Ruminococcus 1` = Species_WGS$s__Ruminococcus_bromii[match(rownames(Genus_WGS), rownames(Species_WGS))]

# subset WGS
Genus_WGS = as.data.frame(t(Genus_WGS))
Genus_WGS = Genus_WGS[which(rownames(Genus_WGS) %in% tumor.RA$genus),]

# subset 16s to same genera in WGS
tumor.RA = tumor.RA[which(tumor.RA$genus %in% rownames(Genus_WGS)),]
rownames(tumor.RA) = tumor.RA$genus
tumor.RA$genus = NULL

# subset 16s to same patients as WGS
tumor.RA = tumor.RA[,colnames(Genus_WGS)]

# convert to matrix 
tumor.RA = as.matrix(t(tumor.RA))
Genus_WGS = as.matrix(t(Genus_WGS))

# reorder the same 
colnames(tumor.RA) == colnames(Genus_WGS)
Genus_WGS = Genus_WGS[,colnames(tumor.RA)]
colnames(tumor.RA) == colnames(Genus_WGS)

# correlation 
matrix.cor = cor(tumor.RA, Genus_WGS, method = "spearman")

# extract diag
diag.correlation = as.data.frame(diag(matrix.cor))

median(diag.correlation$`diag(matrix.cor)`)
mean(diag.correlation$`diag(matrix.cor)`)

# write.csv(diag.correlation, file = "./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/WGS_validation/spearman_correlation_16s_tumor_WGS_common_genera_from_model.csv")

# df with model genera
df = data.frame(genera = Best_CV_coef$covariates, correlation = NA)

# add correlation to the availabe 
df$correlation = diag.correlation$`diag(matrix.cor)`[match(df$genera, rownames(diag.correlation))]
rownames(df) = df$genera
df$genera = NULL

# diag.correlation = as.matrix(t(diag.correlation))

df = as.matrix(t(df))
#rownames(diag.correlation)[1] = "correlation"

#colnames(diag.correlation)

#diag.correlation = data.frame(diag.correlation)
#colnames(diag.correlation) = gsub("\\.", " ", colnames(diag.correlation))
#colnames(diag.correlation)[-which(colnames(diag.correlation) %in% Best_CV_coef$covariates)]
df = df[, Best_CV_coef$covariates]
df = as.data.frame(t(df))
df = as.matrix(df)
rownames(df)[1] = "correlation"

colnames(df) == Best_CV_coef$covariates

# heatmap colors
#col_fun = colorRamp2(c(min(diag.correlation), 0, max(diag.correlation)), c("#FDF5E6", "#FFDAB9", "red"))
col_fun = colorRamp2(c(min(diag.correlation), 0, max(diag.correlation)), c("#00BFFF", "#87CEFA", "#00008B"))

# plot heatmap
HM = Heatmap(df,
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 8),
             cluster_rows = F, cluster_columns = F,
             column_dend_reorder = F,
             row_dend_reorder = F,
             #left_annotation = hr,
             show_column_names = T, 
             #top_annotation = col.an, 
             col = col_fun,
             na_col = "#ECECEC",
             show_heatmap_legend = TRUE,
             row_names_max_width = unit(6, "in"))

svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_16s_tumor_WGS_41_v2.svg"),
    width = 10, height = 2.1)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_16s_tumor_WGS_41_v2.pdf"),
    width = 10, height = 2.1)

HM = draw(HM, heatmap_legend_side = "left")
dev.off()

