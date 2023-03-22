
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "easyGgplot2")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
geneset = "ESTIMATE"

# Load data
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
load("./Analysis/008_ESTIMATE/Biolinks_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")

# Analysis
table_cluster_assignment$Patient_ID = substring(rownames(table_cluster_assignment), 1, 12)
table_cluster_assignment$CMS = Rfcms$RF.predictedCMS[match(table_cluster_assignment$Patient_ID, substring(rownames(Rfcms), 1, 12))]
table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(filtered.norm.RNAseqData)),]

plot_df = data.frame(Patient_ID = table_cluster_assignment$Patient_ID, ES = NA, ICR_cluster = table_cluster_assignment$ICR_HML, 
                     CMS = table_cluster_assignment$CMS)
plot_df$ICR_cluster = factor(plot_df$ICR_cluster, levels =c("ICR Low", "ICR Medium", "ICR High"))
plot_df$CMS = factor(plot_df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))

dir.create("./Figures/024_violin_ESTIMATE", showWarnings = FALSE)

plot_df = plot_df[-which(is.na(plot_df$CMS)),]

i=1
for (i in 1:ncol(ESTIMATE)){
  var = colnames(ESTIMATE)[i]
  plot_df$ES = ESTIMATE[,var][match(plot_df$Patient_ID, substring(rownames(ESTIMATE), 1, 12))]
  plot = ggplot(plot_df, aes(x=ICR_cluster, y = ES, fill = ICR_cluster)) +
    facet_grid(.~CMS) +
    geom_rect(aes(fill = CMS), xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3) +
    scale_fill_manual(values = c("ICR High" = alpha("red", 1), 
                                 "ICR Medium" = alpha("green", 1),
                                 "ICR Low" = alpha("blue", 1),
                                 "CMS1" = alpha("#FFD09E", 1),
                                 "CMS2" = alpha("#97BDD8", 1), 
                                 "CMS3" = alpha("#F4C0D5", 1), 
                                 "CMS4" = alpha("#99CFBE", 1))) +
    geom_violin() +
    geom_boxplot(width=.1, outlier.shape = NA) +
    #geom_jitter(width = 0.15, size = 0.2) +
    #theme_bw() +
    #stat_compare_means(method = "t.test", label = "p.signif",
    #comparisons = list(c("ICR High", "ICR Low"),
    #c("ICR Medium", "ICR Low"),
    #c("ICR High", "ICR Medium"))) +
    ylab(paste0(var)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 10),
          strip.background = element_blank(),
          strip.text = element_blank(),
          aspect.ratio = 1.3/1,
          legend.position = "none") 
  
  png(paste0("./Figures/024_violin_ESTIMATE/", var, "_violinplot.png"),
      width = 6, height = 2, units = "in", res = 600)
  plot(plot)
  dev.off()
}
