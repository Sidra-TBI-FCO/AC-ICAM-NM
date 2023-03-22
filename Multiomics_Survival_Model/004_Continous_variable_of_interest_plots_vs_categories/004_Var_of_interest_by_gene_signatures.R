

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ComplexHeatmap", "circlize", "Hmisc"))

# Set parameters
var_of_interest = "var15_Microbiome_risk_score"
subset = "246"

# Load data
load("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

if(subset == "246"){
  df = df[which(df$Microbiome_Paired_TN == "yes"),]
}

# Set parameters
dir.create(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots"), showWarnings = FALSE)

# Prepare data
immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in% 
                                      df$Patient_ID),]

immune_sig_df$variable = df[, var_of_interest][match(substring(rownames(immune_sig_df), 1, 3),
                                                     df$Patient_ID)]

# Calculate correlation
mat_cor = cor(immune_sig_df, method = "spearman")
mat_cor = as.matrix(mat_cor)

corr = rcorr(as.matrix(immune_sig_df),type="spearman")
mat_cor = corr$r
mat_p_value = corr$P

# Heatmap
mat_cor = mat_cor[1:103,104, drop = FALSE]
mat_p_value = mat_p_value[1:103,104, drop = FALSE]

rownames(mat_cor) = gsub(".*\\-", "", rownames(mat_cor))
rownames(mat_p_value) = gsub(".*\\-", "", rownames(mat_p_value))

col_fun = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "#ff2626"))

png(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots/", var_of_interest, "_correlation_with_immune_gene_signatures_", 
           "in_ACICAM", subset, ".png"), 
    res = 600, width = 2.5, height = 7, units = "in")
HM = Heatmap(mat_cor,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5),
             column_title_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
             cluster_rows = TRUE, cluster_columns = T,
             show_column_names = T,
             col = col_fun,
             show_row_dend = FALSE,
             #show_heatmap_legend = FALSE,
             #col = color,
             heatmap_legend_param =list(title_gp=gpar(fontsize=10, fontface="bold"),legend_width=unit(8,"cm"),legend_position = "left"),
             #row_names_max_width = unit(5, "cm")
             # theme(legend.position = "none")
             row_names_max_width = unit(6, "in")
)
HM = draw(HM)

dev.off()

order = row_order(HM)

mat_p_value = mat_p_value[order,]
mat_p_value = mat_p_value < 0.05 

col_fun2 = colorRamp2(c(-0.3, 0, 0.30), c("blue", "white", "darkorange"))

png(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots/p_value_", var_of_interest, "_correlation_with_immune_gene_signatures_", 
           "in_ACICAM", subset, ".png"), 
    res = 600, width = 2.5, height = 7, units = "in")
HM = Heatmap(mat_p_value,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5),
             column_title_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
             cluster_rows = FALSE, cluster_columns = T,
             show_column_names = T,
             col = col_fun2,
             show_row_dend = FALSE,
             #show_heatmap_legend = FALSE,
             #col = color,
             heatmap_legend_param =list(title_gp=gpar(fontsize=10, fontface="bold"),legend_width=unit(8,"cm"),legend_position = "left"),
             #row_names_max_width = unit(5, "cm")
             # theme(legend.position = "none")
             row_names_max_width = unit(6, "in")
)
HM = draw(HM)

dev.off()




