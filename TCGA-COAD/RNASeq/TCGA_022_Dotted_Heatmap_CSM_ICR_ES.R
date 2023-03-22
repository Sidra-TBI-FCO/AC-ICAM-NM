
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize",
                                   "dendsort", "stringr")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
geneset = "ConsensusTME_COAD"

# Manual adjustment
order = c("ICR score", "Fibroblasts", "Endothelial", "Neutrophils", "Monocytes", "Macrophages", "Eosinophils", "Mast cells",
          "Cytotoxic cells", "NK cells", "T regulatory cells", "T cells gamma delta", "T cells CD8", "T cells CD4",
          "Immune Score", "Macrophages M1", "Dendritic cells", "Macrophages M2", "Plasma cells", "B cells")

# Load data
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
load(paste0("./Analysis/Deconvolution_and_GSEA/", geneset,"_ES.Rdata"))
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")

# Analysis
plot_df = data.frame(t(ES))
plot_df$ICR = table_cluster_assignment$ICR_HML[match(rownames(plot_df),
                                                     rownames(table_cluster_assignment))]
plot_df$ICR_score = table_cluster_assignment$ICRscore[match(rownames(plot_df),
                                                            rownames(table_cluster_assignment))]
plot_df$CMS = Rfcms$RF.predictedCMS[match(rownames(plot_df),
                                          rownames(Rfcms))]
plot_df$ICR = factor(plot_df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))
plot_df_agg = aggregate(.~CMS+ICR, plot_df, FUN = median)
rownames(plot_df_agg) = paste(plot_df_agg$ICR, plot_df_agg$CMS)
plot_df_agg$ICR = NULL
plot_df_agg$CMS = NULL
plot_df_agg = as.matrix(t(plot_df_agg))

rownames(plot_df_agg) = gsub("\\_", " ", rownames(plot_df_agg))
plot_df_agg = plot_df_agg[order,]

# z-score Expression.matrix
plot_df_agg_z = plot_df_agg
for(j in 1: nrow(plot_df_agg_z))  {
  plot_df_agg_z[j,] = (plot_df_agg[j,]-mean(plot_df_agg[j,]))/sd(plot_df_agg[j,]) # z-score the enrichment matrix
}

ha = HeatmapAnnotation(CMS = str_sub(colnames(plot_df_agg_z), -4, -1),
                       `ICR cluster` = gsub(" CMS*.", "", colnames(plot_df_agg_z)),
                       col = list(`ICR cluster` = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  `CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74")))
dir.create("./Figures/022_Dotted_Heatmap_ES", showWarnings = FALSE)

rownames(plot_df_agg_z) = gsub("\\.", " ", rownames(plot_df_agg_z))

#dend = dendsort(hclust(dist(plot_df_agg_z)))

col_fun = colorRamp2(c(-2.7,0, 2.4), c("blue","white", "red"))

png(paste0("./Figures/022_Dotted_Heatmap_ES/Dotted_heatmap_", geneset, ".png"),
    res = 600, width = 7.5, height = 4.5, units = "in")
Heatmap(plot_df_agg_z, cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Mean enrichment score (z-scored)",
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = "#EBEBEB"))
          grid.circle(x = x, y = y, r = 0.024,gp = gpar(fill = col_fun(plot_df_agg_z[i, j]), col = NA))
        }
)
dev.off()
