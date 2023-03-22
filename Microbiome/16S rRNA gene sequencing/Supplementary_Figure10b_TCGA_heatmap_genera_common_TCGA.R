
# spearman correlation 
# ordered as AC-ICAM

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))

source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "corrplot", "scatterplot3d", "ggplot2", "RColorBrewer", "stringr", "ComplexHeatmap", "stringr","circlize"))

# load data
load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_filtered_per_10_RA_0.01.Rdata")
df = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/genera_annotation_file_updated.csv")
ac.icam = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/spearman_common_genera/AC_ICAM_hclust_spearman_correlation_ordred_clusters_v1_OCT_filtered_TCGA.csv")
#cor.matrix = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/Spearman_correlation_138x138_genera_in_tumor.csv")
load("./Analysis/Exploration_reviewers/Microbiome/round_two/sparcc_selva_genus_annotated_with_clusters_corrplot_500_irt.Rdata") 

## TCGA ####
# new object for microbiome
microbiome.TCGA = tcga_bacteria_df


df$genus = str_trim(df$genus, "both")
df$phylum = str_trim(df$phylum, "both")
df$family = str_trim(df$family, "both")

microbiome.TCGA = microbiome.TCGA[which(rownames(microbiome.TCGA) %in% df$genus),]

microbiome.TCGA = t(microbiome.TCGA)
microbiome.TCGA = as.matrix(microbiome.TCGA)

# correlation AC-ICAM
TCGA.cor= cor(microbiome.TCGA, method = "spearman")
TCGA.plot = corrplot(TCGA.cor, tl.col = "black", order = "hclust", method = "color", tl.cex = 0.7/par("cex"), cl.cex = 0.7/par("cex"))

# table(df$phylum)
# table(df$family)
# unique(df$family)
# unique(df$phylum)

# plot correlation using heatmap
cor.plot = TCGA.plot$corr
class(cor.plot)

ac.icam = ac.icam[-grep("Ruminococcus 1", ac.icam$x),] 
ac.icam$x[grep("Ruminococcus 2", ac.icam$x)] = "Ruminococcus"

cor.plot = cor.plot[ac.icam$x, ]
cor.plot = cor.plot[,ac.icam$x]

df = df[which(df$genus %in% rownames(cor.plot)),]

rownames(df) = df$genus
df = df[rownames(cor.plot),]

table(df$cluster)

df$cluster = factor(df$cluster, levels = c("cluster1", "cluster2" ,"cluster3", "cluster4" ,"cluster5" , "cluster7"))


df$risk.group = factor(df$risk.group, levels = c("Low" , "High" , "Not included"))

clusters_df = sparcc_corr
rownames(clusters_df) = clusters_df$genus

clusters_df = clusters_df[-grep("Ruminococcus 1", clusters_df$genus),]
clusters_df$genus = gsub("Ruminococcus 2","Ruminococcus" ,clusters_df$genus)


clusters_df = clusters_df[df$taxanomy,]

table(clusters_df$cluster)
clusters_df$cluster = factor(clusters_df$cluster, levels = c("cluster1", "cluster2", "cluster4","cluster5" ,"cluster8", "No cluster"))


# v3
hr = rowAnnotation(df = data.frame(SparCC = clusters_df$cluster, Spearman = df$cluster, MBR.group = df$risk.group),
                   col = list(SparCC = c("cluster1" = "#FF1493", "cluster2" = "#0000FF", "cluster4" = "#48D1CC", 
                                         "cluster5" = "#800000", "cluster8" = "#FFA500", "No cluster" = "grey"),
                              Spearman = c("cluster1" = "#ED1B33", "cluster2" = "#C0DEA1", "cluster3" = "#20B2AA", "cluster4" = "#FFD28A",
                                           "cluster5" = "#CB56A0",
                                           "cluster7" = "#242062"),

                              MBR.group = c("Low" = "#3F55A4", "High" = "#F06476", "Not included" = "grey")),
                   #gp = gpar(col = "#DCDCDC"),
                   show_legend = TRUE)


# plot heatmap
HM = Heatmap(cor.plot,
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
svg(paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/common_genera/heatmap_TCGA_spearman_annotated_spearman_sparcc_clusters_common_TCGA_v1.svg"),
    width = 7, height = 5)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

### pdf 

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/common_genera/heatmap_TCGA_spearman_annotated_spearman_sparcc_clusters_common_TCGA_v1.pdf"),
    width = 7, height = 5)

HM = draw(HM, heatmap_legend_side = "left")
dev.off() 

