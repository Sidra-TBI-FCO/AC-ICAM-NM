

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize",
                                   "dendsort", "stringr")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
geneset = "ConsensusTME_COAD" # "ConsensusTME_COAD", "Bindea_ORIG"
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load(paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/", geneset,"_ES.Rdata"))
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Analysis
ES = ES[, -which(substring(colnames(ES), 1, 3) %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)])]
plot_df = data.frame(t(ES))
plot_df$ICR = table_cluster_assignment$ICR_HML[match(rownames(plot_df),
                                                     rownames(table_cluster_assignment))]
plot_df$ICRscore = table_cluster_assignment$ICRscore[match(rownames(plot_df),
                                                     rownames(table_cluster_assignment))]
plot_df$CMS = Rfcms$RF.predictedCMS[match(rownames(plot_df),
                                          rownames(Rfcms))]
plot_df$ICR = factor(plot_df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))
plot_df$CMS = factor(plot_df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4")) # reverse order because plot will be mirrored
plot_df_agg = aggregate(.~CMS+ICR, plot_df, FUN = median)
rownames(plot_df_agg) = paste(plot_df_agg$ICR, plot_df_agg$CMS)
plot_df_agg$ICR = NULL
plot_df_agg$CMS = NULL
plot_df_agg = as.matrix(t(plot_df_agg))

# z-score Expression.matrix
plot_df_agg_z = plot_df_agg
for(j in 1: nrow(plot_df_agg_z))  {
  plot_df_agg_z[j,] = (plot_df_agg[j,]-mean(plot_df_agg[j,]))/sd(plot_df_agg[j,]) # z-score the enrichment matrix
}

ha = HeatmapAnnotation(CMS = str_sub(colnames(plot_df_agg_z), -4, -1),
                       `ICR cluster` = gsub(" CMS*.", "", colnames(plot_df_agg_z)),
                       col = list(`ICR cluster` = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  `CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74")))
dir.create("./Figures/Trimmed_p/022_Dotted_Heatmap_ES", showWarnings = FALSE)

rownames(plot_df_agg_z) = gsub("\\.", " ", rownames(plot_df_agg_z))
rownames(plot_df_agg_z) = gsub("\\_", " ", rownames(plot_df_agg_z))

#dend = dendsort(hclust(dist(plot_df_agg_z)))

col_fun = colorRamp2(c(-2,0, 2.13), c("blue","white", "red"))

pdf(paste0("./Figures/Trimmed_p/022_Dotted_Heatmap_ES/Fig1c__ICR_continuous_Dotted_heatmap_", geneset, ".pdf"),
    width = 8.0, height = 4.5) 
Heatmap(plot_df_agg_z, cluster_rows = TRUE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Mean enrichment score (z-scored)",
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = "#EBEBEB"))
          grid.circle(x = x, y = y, r = 0.025,gp = gpar(fill = col_fun(plot_df_agg_z[i, j]), col = NA))
        }
)
dev.off()
