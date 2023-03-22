
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
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

TCR_Overview$HLM_cluster = factor(TCR_Overview$HLM_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
levels(TCR_Overview$HLM_cluster) = c("Low", "Medium", "High")
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
TCR_Overview$CMS = Rfcms$RF.predictedCMS[match(TCR_Overview$Patient_ID, substring(rownames(Rfcms), 1, 3))]
TCR_Overview$CMS = factor(TCR_Overview$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
TCR_Overview = TCR_Overview[-which(is.na(TCR_Overview$CMS)),]

# Plotting
df_plot = TCR_Overview[, c("ICRscore", "CMS", variable, "HLM_cluster")]

df_plot$ICRcluster = df_plot$HLM_cluster
df_plot$HLM_cluster = NULL

plot = ggplot(df_plot, aes(x = ICRcluster, y = productive_clonality)) + 
  facet_grid(~CMS) +
  geom_boxplot(width=1, outlier.shape = NA, aes(fill = ICRcluster)) +
  geom_jitter(size = 0.8, width = 0.1) +
  scale_fill_manual(values = c("High" = "red", "Medium" = "green",
                               "Low" = "blue")) +
  xlab("") +
  ylab("TCR productive clonality \nImmunoSeq") +
  theme_bw() +
  ylim(0, 0.4) +
  #stat_compare_means(method = "t.test", comparisons = list(c("Medium", "Low"),
   #                                                        c("High", "Medium"),
    #                                                       c("High", "Low")),
     #                label = "p.signif") +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12))

pdf(file = "./Figures/TCR/013_Boxplot_clonality_CMS/EDF4a_Violin_plot_clonality_by_ICR_facet_by_CMS_stats.pdf",
    width = 6, height = 3)
plot(plot)
dev.off()
