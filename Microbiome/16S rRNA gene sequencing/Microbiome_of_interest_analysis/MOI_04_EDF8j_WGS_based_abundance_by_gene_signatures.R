

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ComplexHeatmap", "circlize", "Hmisc"))

# Set parameters
Species = "s__Fusobacterium_nucleatum"
All_signatures = "" # "All_signatures" or ""
MSI = "MSI-H" # or "MSI-H"

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

if(All_signatures == "All_signatures"){
  load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_All_Immune_gene_signatures_table.Rdata")
}

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures"), showWarnings = FALSE)

# Getting MSS and MSI-H
include_patients = MANTIS$Patient_ID[which(MANTIS$MSI == MSI)]

# Prepare data
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in% 
                                      substring(colnames(Species_WGS), 1, 3)),]

immune_sig_df$Abundance = Species_WGS[Species,][match(substring(rownames(immune_sig_df), 1, 3),
                                                             substring(colnames(Species_WGS), 1, 3))]

save(immune_sig_df, file = paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/004_correlation_with_gene_signatures/immune_sig_df_",
                                  Species, ".Rdata"))

if(MSI == ""){}else{
  immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in% include_patients),]
}

# Calculate correlation
mat_cor = cor(immune_sig_df, method = "spearman")
mat_cor = as.matrix(mat_cor)

corr = rcorr(as.matrix(immune_sig_df),type="spearman")
mat_cor = corr$r
mat_p_value = corr$P

# Heatmap
mat_cor = mat_cor[1:(ncol(mat_cor)-1), ncol(mat_cor), drop = FALSE]
mat_p_value = mat_p_value[1:(ncol(mat_p_value)-1), ncol(mat_p_value), drop = FALSE]

mat_FDR_value = mat_p_value
mat_FDR_value  = data.frame(mat_FDR_value)
mat_FDR_value$FDR = p.adjust(mat_FDR_value$Abundance, method = "BH", nrow(mat_FDR_value))

test = cbind(mat_cor, mat_FDR_value)
colnames(test) = c("Rho", "p_value", "FDR_BH")

write.csv(test, file = paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/", 
                              "004_correlation_with_gene_signatures/Stats_", MSI, "_", Species, "_WGS_correlation_EDF8.csv"),
          row.names = TRUE)

# Generate overview with NK cells
#NK_sub = mat_cor[grep("NK", rownames(mat_cor)), 1, drop = FALSE]
#NK_sub_p_value = mat_p_value[grep("NK", rownames(mat_p_value)), 1, drop = FALSE]
#NK_sub = data.frame(NK_sub)
#NK_sub$p_value = NK_sub_p_value
#colnames(NK_sub) = c("Spearman_correlation", "p_value")

#dir.create("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/004_correlation_with_gene_signatures", showWarnings = FALSE)
#write.csv(NK_sub, file = paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/004_correlation_with_gene_signatures/",
#                               Genus, "_by_NK_signatures.csv"))

rownames(mat_cor) = gsub(".*\\-", "", rownames(mat_cor))
rownames(mat_p_value) = gsub(".*\\-", "", rownames(mat_p_value))

col_fun = colorRamp2(c(-0.3, 0, 0.30), c("blue", "white", "#ff2626"))

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures/Supplementary_Figure11e_WGS_based_", Species, "_", MSI,"_", All_signatures, "_correlation_with_immune_gene_signatures.pdf"), 
    width = 2.5, height = 7)
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures/Supplementary_Figure11e_WGS_based_", Species, "_", MSI,"_", All_signatures, "_correlation_pval_with_immune_gene_signatures.pdf"), 
    width = 2.5, height = 7)
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




