# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "corrplot", "ggplot2", "RColorBrewer", "stringr", "ComplexHeatmap","circlize", "openxlsx"))

# load data
load("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/matrices/sparcc_selva_hclust_matrix.Rdata")
annotation = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/genera_annotation_file_updated.csv")
load("./Analysis/Exploration_reviewers/Microbiome/round_two/sparcc_selva_genus_annotated_with_clusters_corrplot_500_irt.Rdata") 
load("./Analysis/Exploration_reviewers/Microbiome/round_two/old/SparCC/genera_common_AC_ICAM_filtered_TCGA.Rdata")

df = df[order(df$cluster),]

# set rownames of cor.matrxi 
cor.matrix = sparcc.matrix

# order cor.matrix
cor.matrix = cor.matrix[df$taxanomy,]
cor.matrix = cor.matrix[,df$taxanomy]

# set as matrix
cor.matrix = as.matrix(cor.matrix)


table(df$cluster)
df$cluster = factor(df$cluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster7"))

# order clusters
clusters_df = sparcc_corr
clusters_df = clusters_df[which(clusters_df$genus %in% df$taxanomy),]

rownames(clusters_df) = clusters_df$genus
clusters_df = clusters_df[rownames(cor.matrix),]

table(clusters_df$cluster)
clusters_df$cluster = factor(clusters_df$cluster, levels = c("cluster1", "cluster2", "cluster4","cluster5" ,"cluster8", "No cluster"))

table(df$risk.group)
df$risk.group = factor(df$risk.group, levels = c("Low", "High", "Not included"))


rownames(cor.matrix) = gsub(".*\\D_5__", "", rownames(cor.matrix))


# v3
hr = rowAnnotation(df = data.frame(SparCC = clusters_df$cluster, Spearman = df$cluster, MBR.group = df$risk.group),
                   col = list(SparCC = c("cluster1" = "#FF1493", "cluster2" = "#0000FF", "cluster4" = "#48D1CC", 
                                         "cluster5" = "#800000", "cluster8" = "#FFA500", "No cluster" = "grey"),
                              Spearman = c("cluster1" = "#ED1B33", "cluster2" = "#C0DEA1", "cluster3" = "#20B2AA", "cluster4" = "#FFD28A",
                                           "cluster5" = "#CB56A0", "cluster7" = "#242062"),

                              MBR.group = c("Low" = "#3F55A4", "High" = "#F06476", "Not included" = "grey")),
                   #gp = gpar(col = "#DCDCDC"),
                   show_legend = TRUE)


# plot heatmap
HM = Heatmap(cor.matrix,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 10),
             column_title_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6.5),
             cluster_rows = F, cluster_columns = F, 
             left_annotation = hr,
             column_dend_reorder = F,
             row_dend_reorder = F,
             show_column_names = F, 
             show_row_names = T,
             #col = col_fun,
             show_heatmap_legend = TRUE,
             na_col="grey",
             row_names_max_width = unit(6, "in")
) 

#dir.create("./Figures/Microbiome/017_Spearman_Genus_T_Immune_sign_Heatmap", showWarnings = FALSE)
svg(paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/common_genera/heatmap_sparcc_annotated_spearman_sparcc_MBR_clusters_common_TCGA_v1.svg"),
    width = 7, height = 5)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

### pdf 

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/common_genera/heatmap_sparcc_annotated_spearman_sparcc_MBR_clusters_common_TCGA_v1.pdf"),
    width = 7, height = 5)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

