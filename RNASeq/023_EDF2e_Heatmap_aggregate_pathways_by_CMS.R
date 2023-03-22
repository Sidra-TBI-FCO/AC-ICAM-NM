
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize",
                                   "dendsort")                                                                   
ibiopak(required.bioconductor.packages)
ipak("stringr")

# Set parameters
geneset = "Selected.Pathways"
Tolga_pathway_excluded = "Tolga_pathway_excluded"
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/", geneset,"_ES.Rdata"))
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Analysis
ES = ES[, -which(substring(colnames(ES), 1, 3) %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)])]
plot_df = data.frame(t(ES))

plot_df$`TGFB PCA 17349583` = immune_sig_df$`Expression Signature - TGFB PCA 17349583`[match(rownames(plot_df),
                                                                                      rownames(immune_sig_df))]
plot_df$`TGFB score 21050467` = immune_sig_df$`Expression Signature - TGFB score 21050467`[match(rownames(plot_df),
                                                                                               rownames(immune_sig_df))]
plot_df$`Module11 Prolif score` = immune_sig_df$`Expression Signature - Module11 Prolif score`[match(rownames(plot_df),
                                                                                                 rownames(immune_sig_df))]
plot_df$`Angiogenesis` = immune_sig_df$`Expression Signature - Angiogenesis`[match(rownames(plot_df),
                                                                                                     rownames(immune_sig_df))]
plot_df$X.HM..TGF.beta.signaling = NULL
plot_df$X.LM..Proliferation = NULL
plot_df$X.HM..Angiogenesis = NULL

plot_df$ICR = table_cluster_assignment$ICR_HML[match(rownames(plot_df),
                                                     rownames(table_cluster_assignment))]
plot_df$CMS = Rfcms$RF.predictedCMS[match(rownames(plot_df),
                                          rownames(Rfcms))]
plot_df$ICR = factor(plot_df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))
plot_df_agg = aggregate(.~CMS, plot_df, FUN = median)
rownames(plot_df_agg) = plot_df_agg$CMS
plot_df_agg$ICR = NULL
plot_df_agg$CMS = NULL
plot_df_agg = as.matrix(t(plot_df_agg))

# z-score Expression.matrix
plot_df_agg_z = plot_df_agg
for(j in 1: nrow(plot_df_agg_z))  {
  plot_df_agg_z[j,] = (plot_df_agg[j,]-mean(plot_df_agg[j,]))/sd(plot_df_agg[j,]) # z-score the enrichment matrix
}

plot_df_agg_z = plot_df_agg_z[-grep("Breast", rownames(plot_df_agg_z)),]

if(Tolga_pathway_excluded == "Tolga_pathway_excluded"){
  plot_df_agg_z = plot_df_agg_z[-grep("X.TPW..", rownames(plot_df_agg_z)),]
}


ha = HeatmapAnnotation(CMS = colnames(plot_df_agg_z),
                       col = list(`CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74")),
                       show_legend = FALSE,
                       annotation_name_gp = gpar(fontsize = 17))
dir.create("./Figures/Trimmed_p/023_aggr_Heatmap_ES", showWarnings = FALSE)

rownames(plot_df_agg_z) = gsub("X.HM..", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("X.IPA..", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("X.TBI..", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("X.TPW..", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("X.LM..", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("\\.", " ", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("signaling", "sign.", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("Signaling", "sign.", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("Reactive oxigen species", "ROS", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("Epithelial mesenchymal transition", "EMT", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub(" by Telomerase", "", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub(" in Eukaryotes", "", rownames(plot_df_agg_z))
#dend = dendsort(hclust(dist(plot_df_agg_z)))

#col_fun = colorRamp2(c(-1.46,0, 1.48), c("blue","white", "red"))

ha = HeatmapAnnotation(CMS = colnames(plot_df_agg_z),
                       col = list(`CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74")),
                       show_legend = FALSE,
                       annotation_name_gp = gpar(fontsize = 17))
pdf(paste0("./Figures/Trimmed_p/023_aggr_Heatmap_ES/Fig2e_Aggregated_heatmap_", geneset, ".pdf"),
    width = 5.3, height = 11)
Heatmap(plot_df_agg_z, cluster_rows = TRUE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Mean enrichment score (z-scored)",
        show_heatmap_legend = FALSE,
        row_names_max_width = unit(7, "in"),
        rect_gp = gpar(col = "white", lwd = 2),
        row_names_gp = gpar(fontsize = 17)
)
dev.off()
