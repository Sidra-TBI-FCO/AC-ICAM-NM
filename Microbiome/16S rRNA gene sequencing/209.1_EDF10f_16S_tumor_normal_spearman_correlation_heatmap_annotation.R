
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "stringr", "corrplot", "openxlsx", "circlize", "ComplexHeatmap"))

# Load data
#load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
tumor.RA = Genus_full_abundance
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_N_246_samples.Rdata")
normal.RA = Genus_full_abundance
rm(Genus_full_abundance)
#load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
Best_CV_coef = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv", stringsAsFactors = FALSE)

# 
Best_CV_coef$covariates = gsub(".*\\D_5__", "", Best_CV_coef$covariates)
#Best_CV_coef = Best_CV_coef[-which(Best_CV_coef$covariates %in% c("Treponema 2", "Candidatus Soleaferrea", "Mogibacterium", "Caproiciproducens", "Halomonas")),]

normal.RA = normal.RA[which(rownames(normal.RA) %in% Best_CV_coef$coefficients),]
tumor.RA = tumor.RA[which(rownames(tumor.RA) %in% Best_CV_coef$coefficients),]

# rename rows
#rownames(normal.RA)[7] = "uncultured (Prevotellaceae)"
#rownames(tumor.RA)[7] = "uncultured (Prevotellaceae)"

rownames(normal.RA) = gsub(".*\\D_5__", "", rownames(normal.RA))
rownames(tumor.RA) = gsub(".*\\D_5__", "", rownames(tumor.RA))

##
normal.RA = as.matrix(t(normal.RA))
tumor.RA = as.matrix(t(tumor.RA))

matrix.cor = cor(normal.RA, tumor.RA, method = "spearman")

diag.correlation = as.data.frame(diag(matrix.cor))
diag.correlation = as.matrix(t(diag.correlation))

rownames(diag.correlation)[1] = "correlation"

colnames(diag.correlation)

#diag.correlation = data.frame(diag.correlation)
#colnames(diag.correlation) = gsub("\\.", " ", colnames(diag.correlation))
#colnames(diag.correlation)[-which(colnames(diag.correlation) %in% Best_CV_coef$covariates)]
diag.correlation = diag.correlation[, Best_CV_coef$covariates]

diag.correlation = as.data.frame(t(diag.correlation))

colnames(diag.correlation) == Best_CV_coef$covariates

# heatmap colors
col_fun = colorRamp2(c(min(diag.correlation), 0.4, max(diag.correlation)), c("#FDF5E6", "#FFDAB9", "red"))
#col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "#8B008B", "#ff2626"))

# plot heatmap
HM = Heatmap(diag.correlation,
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 8),
             cluster_rows = F, cluster_columns = F,
             column_dend_reorder = F,
             row_dend_reorder = F,
             #left_annotation = hr,
             show_column_names = T, 
             #top_annotation = col.an, 
             col = col_fun,
             show_heatmap_legend = TRUE,
             row_names_max_width = unit(6, "in"))

svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_normal_tumor_41_v1.svg"),
    width = 10, height = 2.1)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 


pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_normal_tumor_41_v1.pdf"),
    width = 10, height = 2.1)

HM = draw(HM, heatmap_legend_side = "left")
dev.off()
