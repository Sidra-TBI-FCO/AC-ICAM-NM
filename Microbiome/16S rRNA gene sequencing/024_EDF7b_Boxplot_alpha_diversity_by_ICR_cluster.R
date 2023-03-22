
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "stringr"))

# Set parameters
#Tissue = "T"

# Load data
load("./Processed_Data/Microbiome/23_Alpha_diversity/Sample_based_filted/alpha_diversity_TN_pairs.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

# Only select Tumor or normal samples
#df_alpha = df_alpha[which(df_alpha$Tissue == Tissue),]
df_alpha$ICR_cluster = table_cluster_assignment$ICR_HML[match(substring(df_alpha$Sample_ID, 1, 3),
                                                              substring(rownames(table_cluster_assignment), 1, 3))]
df_alpha$ICR_cluster = factor(df_alpha$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

vars = colnames(df_alpha)[4:7]

i=1
for (i in 1:length(vars)){
  var = vars[i]
  plot = ggplot(df_alpha, aes(x = ICR_cluster, y = get(var), fill = ICR_cluster)) +
    scale_fill_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
    facet_grid(.~Tissue) +
    geom_boxplot(outlier.shape = NA) +
    #ggtitle(paste0(Tissue)) +
    geom_jitter(size = 0.5, width = 0.1) +
    theme_bw() +
    xlab("") +
    ylab(var) +
    stat_compare_means(method = "wilcox", comparisons = list(c("ICR High", "ICR Medium"),
                                                             c("ICR Medium", "ICR Low"),
                                                             c("ICR High", "ICR Low"))) +
    theme(axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_text(size = 15, color = "black"),
          #axis.text.x = element_blank(),
          axis.text.x = element_text(size = 15, color = "black",
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 15, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(size = 15, color = "black"),
          legend.position = "none")
  
  dir.create("./Figures/Microbiome/024_Boxplot_alpha_diversity_by_ICR_cluster", showWarnings = FALSE)
  
  #png(paste0("./Figures/Microbiome/024_Boxplot_alpha_diversity_by_ICR_cluster/Facet_Boxplot_", var, "_by_ICR_tissue.png"),
   #   res = 600, width = 4, height = 4, units = "in")
  #plot(plot)
  #dev.off()
  
  pdf(paste0("./Figures/Microbiome/024_Boxplot_alpha_diversity_by_ICR_cluster/Facet_Boxplot_", var, "_by_ICR_tissue.pdf"),
      width = 4, height = 5)
  plot(plot)
  dev.off()
}

