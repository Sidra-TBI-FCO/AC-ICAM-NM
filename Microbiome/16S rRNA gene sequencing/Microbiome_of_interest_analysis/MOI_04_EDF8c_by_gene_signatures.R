
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ComplexHeatmap", "circlize", "Hmisc"))

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
# rownames(Genus_full_abundance)[grep("Fuso", rownames(Genus_full_abundance))
Genus = "D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
All_signatures = "" # "All_signatures" or ""

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

if(All_signatures == "All_signatures"){
  load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_All_Immune_gene_signatures_table.Rdata")
}

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures"), showWarnings = FALSE)

# Prepare data
immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in% 
                                      substring(colnames(Genus_full_abundance), 1, 3)),]

immune_sig_df$Abundance = Genus_full_abundance[Genus,][match(substring(rownames(immune_sig_df), 1, 3),
                                                             substring(colnames(Genus_full_abundance), 1, 3))]

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
                              "004_correlation_with_gene_signatures/Stats_", Genus, "_16S_correlation_EDF8.csv"),
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures/EDF8c_", Genus, All_signatures, "_correlation_with_immune_gene_signatures.pdf"), 
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_04_Genus_correlation_all_signatures/EDF8c_", Genus, All_signatures, "_correlation_pval_with_immune_gene_signatures.pdf"), 
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

