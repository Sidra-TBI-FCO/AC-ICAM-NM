
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "reshape2", "ComplexHeatmap", "circlize", "dendsort"))

# Set parameters
correlation = "" #"" # significant

# load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_model_genera_risk_score_immune_sig.Rdata") # spearman 
#load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/pearson_correlation/pearson_model_genera_risk_score_immune_sig.Rdata") #pearson
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/HR_tables/013_Aug_2022_OS_HR_table_immune_signatures__cutoff_20_v2.Rdata")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/HR_tables/013_Aug_2022_PFS_HR_table_immune_signatures__cutoff_20_v2.Rdata")
#genera.annotation = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/genus_full_abundance_138_annotated.csv")
load("./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/Immune_signatures_clusters_Roelands.Rdata")
# 42 MBR model
#MBR_CV_coef = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv", stringsAsFactors = FALSE)
#MBR_CV_coef = MBR_CV_coef[,-1]
#colnames(MBR_CV_coef) = c("covariates","coefficients","colors")  
#rownames(MBR_CV_coef) = gsub(".*D_5__", "", MBR_CV_coef$covariates)
#save(MBR_CV_coef,file = "./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.RData")
load("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.RData")


# fix signatures names 
results_all$signature = gsub("Expression Signature - ", "", results_all$signature)
results_all$signature = gsub("Attractor Metagene - ", "", results_all$signature)
results_all$signature = gsub("ConsensusTME - ", "", results_all$signature)
results_all$Name[which(results_all$Name == "risk.score")] = "Risk score"

# new object for results all
annotated = results_all
annotated = annotated[!is.na(annotated$p_val),]
annotated = annotated[!is.na(annotated$rho),]

# genera annotation annotated
annotated = annotated[which(annotated$signature %in% os.hr$Signature),]
unique(annotated$signature)

#annotated = annotated[which(annotated$Name %in% c("Risk score", "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2")),]

#drop Risk score
annotated = annotated[-which(annotated$Name %in% c("Risk score")),]

# calculate FDR
annotated$FDR = NA
annotated$FDR = p.adjust(annotated$p_val, method = "fdr", n = nrow(annotated))

#blank the FDR
#$rho[annotated$Name!="D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2" & annotated$FDR>0.1] = 0
#annotated$rho[annotated$FDR>0.1] = 0
# create correlation matrix
colnames(annotated)[1] = "Genus"
matrix = dcast(annotated, Genus~signature, value.var="rho")
rownames(matrix) = matrix$Genus
matrix$Genus = NULL

class(matrix)
matrix = as.matrix(matrix)

dim(matrix)

# reorder matrix
matrix = matrix[,os.hr$Signature]
matrix = matrix[rev(MBR_CV_coef$covariates),]


# signatures annotation
rownames(annotation) = annotation$Names
annotation = annotation[os.hr$Signature,]

# Set factors
table(annotation$Immune_module)
annotation$Immune_module = factor(annotation$Immune_module, levels = c("Lymphocyte Infiltration", "Macrophage/Monocyte", "IFN Response", "TGF-b Response", "T-cell/Cytotoxic",
                                                                       "Wound Healing", "Unassigned"))

table(annotation$Immune_trait_category)
annotation$Immune_trait_category = factor(annotation$Immune_trait_category, levels = c("ConsensusTME", "Expression Signature", "Attractor Metagene"))

annotation$Module_Roelands[is.na(annotation$Module_Roelands)] = "M7"

table(annotation$Module_Roelands)
annotation$Module_Roelands = factor(annotation$Module_Roelands, levels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7"))

# HR tables
os.hr$pval.cutoff = NA
os.hr$pval.cutoff[which(os.hr$p_value < 0.05)] = "Significant"
os.hr$pval.cutoff[which(os.hr$p_value > 0.05)] = "Not significant"

pfs.hr$pval.cutoff = NA
pfs.hr$pval.cutoff[which(pfs.hr$p_value < 0.05)] = "Significant"
pfs.hr$pval.cutoff[which(pfs.hr$p_value > 0.05)] = "Not significant"

# HR color
hr.col = colorRamp2(c(2.5, 1, 0), c("#FF6666", "white", "#3CB371"))

# heatmap colors
col_fun = colorRamp2(c(min(matrix), 0, max(matrix)), c("blue", "white", "#ff2626"))

# plot legends
# lgd = Legend(col_fun = col_fun, title = "")
# hr_lgd = Legend(col_fun = hr.col, title = "")

# pdf("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/col_fun_legend.pdf", width = 4, height = 4)
# draw(lgd)
# dev.off()
# 
# pdf("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/hr_legend.pdf", width = 4, height = 4)
# draw(hr_lgd)
# dev.off()

rownames(matrix) = gsub(".*\\D_5__", "", rownames(matrix))

#fix genera naming
rownames(matrix)[rownames(matrix) == "uncultured"] = "uncultured Prevotellaceae"
rownames(MBR_CV_coef)[rownames(MBR_CV_coef) == "uncultured"] = "uncultured Prevotellaceae"

# column annotation
col.an = HeatmapAnnotation(df = data.frame(Module = annotation$Module_Roelands, Immune_category = annotation$Immune_trait_category, Immune_module = annotation$Immune_module, HR_OS = os.hr$HR, HR_PFS = pfs.hr$HR),
                           col = list(Module = c("M1" = "#FF1493", "M2" = "#0000FF", "M3" = "#1E90FF", "M4" = "#008080", 
                                                 "M5" = "#800000", "M6" = "#FF8C00", "M7" = "#C0C0C0"),
                                      Immune_category = c("ConsensusTME" = "#C71585", "Expression Signature" = "#2E8B57", "Attractor Metagene" = "#4169E1"),
                                      Immune_module = c("Lymphocyte Infiltration" = "#87CEEB", "Macrophage/Monocyte" = "#90EE90", "IFN Response" = "#FFFF00", 
                                                        "TGF-b Response" = "#9370DB", "T-cell/Cytotoxic" = "#F4A460","Wound Healing" = "#FF7F50", "Unassigned" = "grey"),
                                      HR_OS = hr.col,
                                      HR_PFS = hr.col
                           ), 
                           gp = gpar(col = "#DCDCDC"),
                           show_legend = TRUE)
#row annotation
row.an = rowAnnotation(df = data.frame(direction = MBR_CV_coef$colors),
                       col = list(direction = c("red" = "#3F55A4", "blue" = "#F06476")))
#labels = c("Bad","Good"))

# plot heatmap
HM = Heatmap(matrix,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 10),
             column_title_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
             cluster_rows = F,
             cluster_columns = F,
             #column_split = 4,
             row_split = MBR_CV_coef$colors,
             #left_annotation = hr,
             show_column_names = T, 
             top_annotation = col.an,
             right_annotation = row.an,
             col = col_fun,
             column_dend_reorder = F,
             row_dend_reorder = T,
             show_heatmap_legend = TRUE,
             row_names_max_width = unit(6, "in"))

svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_risk_Score_model_genera_revised_script_v3.svg"),
    width = 16, height = 10)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_correlation_risk_Score_model_genera_revised_script_v3.pdf"),
    width = 16, height = 10)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

# if (correlation == "significant") {
#   
#   annotated$rho[which(annotated$p_val > 0.05)] = NA
#   
#   # create correlation matrix
#   colnames(annotated)[1] = "Genus"
#   matrix = dcast(annotated, Genus~signature, value.var="rho")
#   rownames(matrix) = matrix$Genus
#   matrix$Genus = NULL
#   
#   matrix = as.matrix(matrix)
#   
#   matrix = matrix[,os.hr$Signature]
#   rownames(matrix)[1] = "Ruminococcus 2"
#   
#   rownames()
#   
#   HM = Heatmap(matrix,
#                row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 10),
#                column_title_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
#                cluster_rows = F, cluster_columns = F,
#                #column_split = 4,
#                #row_split = 6,
#                column_dend_reorder = F,
#                row_dend_reorder = F,
#                #left_annotation = hr,
#                show_column_names = T, 
#                #top_annotation = col.an, 
#                col = col_fun,
#                show_heatmap_legend = TRUE,
#                row_names_max_width = unit(6, "in"))
#   
#   svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/04.0_spearman_correlation/spearman_significant_correlation_risk_Score_ruminococcus2_revised_script_v1.svg"),
#       width = 16, height = 10)
#   
#   HM = draw(HM, heatmap_legend_side = "l

#write.csv(MBR_CV_coef,file = "./MBR_genera.csv")
