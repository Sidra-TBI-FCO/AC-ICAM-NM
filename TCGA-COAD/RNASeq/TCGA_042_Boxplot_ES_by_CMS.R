
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "dplyr")
ipak(required.packages)

# Set parameters
Gene.set = "DDR"

# Load data
load(paste0("./Analysis/Deconvolution_and_GSEA/", Gene.set, "_ES.Rdata"))
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")

# Prepare directories
dir.create("./Figures/042_Boxplot_ES_by_CMS", showWarnings = FALSE)
dir.create(paste0("./Figures/042_Boxplot_ES_by_CMS/", Gene.set), showWarnings = FALSE)

# Plot by CMS
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
plot_df = data.frame(Sample = rownames(Rfcms), CMS = Rfcms$RF.predictedCMS, ICR = NA, ES = NA)
plot_df$ICR = table_cluster_assignment$ICR_HML[match(plot_df$Sample, rownames(table_cluster_assignment))]
plot_df$ICR = as.character(plot_df$ICR)
plot_df$ICR[which(plot_df$ICR %in% c("ICR Low", "ICR Medium"))] = "ICR MedLow"
plot_df$ICR = factor(plot_df$ICR, levels = c("ICR MedLow", "ICR High"))

i = 1
for (i in 1:nrow(ES)){
  signature = rownames(ES)[i]
  plot_df$ES = ES[signature,][match(plot_df$Sample, colnames(ES))]
  
  plot = ggplot(plot_df, aes(x = ICR, y = ES)) +
    geom_boxplot(outlier.shape = NA, aes(fill = ICR)) +
    geom_jitter(width = 0.2, size = 0.8) +
    theme_bw() +
    facet_grid(.~CMS, margins = TRUE) +
    ylab(paste0(signature, " \n enrichment score")) +
    theme(axis.text.x = element_text(colour = "black", size = 14, angle = 90, vjust = 1, hjust = 1), 
          axis.text.y = element_text(colour = "black", size = 14),
          axis.title.x = element_text(colour = "black", size = 14),
          axis.title.y = element_text(colour = "black", size = 14),
          strip.background = element_blank()) +
    scale_fill_manual(values = c("ICR MedLow" = "#10B347", "ICR High" = "red"))
  
  png(paste0("./Figures/042_Boxplot_ES_by_CMS/", Gene.set, "/", signature, "_by_CMS_boxplot.png"), 
      res = 600, units = "in", width = 8, height = 4)
  plot(plot)
  dev.off()
}
