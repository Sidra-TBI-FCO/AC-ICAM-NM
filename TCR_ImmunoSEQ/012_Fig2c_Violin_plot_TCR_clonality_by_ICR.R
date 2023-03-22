
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable = "productive_clonality"

# Load data
load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Bindea_ORIG_ES.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")

TCR_Overview$HLM_cluster = factor(TCR_Overview$HLM_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
levels(TCR_Overview$HLM_cluster) = c("Low", "Medium", "High")
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]

# Plotting
df_plot = TCR_Overview[, c("ICRscore", "HLM_cluster", variable)]

plot = ggplot(df_plot, aes(x = HLM_cluster, y = get(variable))) +
  #geom_boxplot(outlier.shape = NA, aes(fill = df_plot$HLM_cluster)) +
  scale_fill_manual(values = c("Low" = "blue", "Medium" = "green","High" = "red")) +
  #geom_jitter(width = 0.2, size = 1.1) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = list(c("Low", "Medium"),
                                                          c("Medium", "High"),
                                                         c("Low", "High"))) +
                   # label = "p.signif") +
  ylab(gsub("_", " ", variable)) +
  xlab("ICR cluster") +
  ylim(0, 0.61)

violin_plot = plot + geom_violin(aes(fill = df_plot$HLM_cluster)) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = df_plot$HLM_cluster)) +
  ylab("") +
  theme(axis.text.x = element_blank()) + xlab("")

pdf(file = "./Figures/TCR/012_Boxplot_clonality_ICR_clusters/Fig2c_Violin_plot_clonality_ICR_cluster_no_stats.pdf",
    width = 3, height = 3.8)
plot(violin_plot)
dev.off()
