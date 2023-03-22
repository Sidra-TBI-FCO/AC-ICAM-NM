
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

plot = ggplot(df_plot, aes(x = CMS, y = productive_clonality)) +
  #geom_boxplot(outlier.shape = NA, aes(fill = df_plot$HLM_cluster)) +
  scale_fill_manual(values = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
                               "CMS3" = "#E97AA8", "CMS4" = "#009E74")) +
  #geom_jitter(width = 0.2, size = 1.1) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.position = "none") +
  stat_compare_means(method = "t.test", comparisons = list(c("CMS1", "CMS2"),
                                                           c("CMS1", "CMS3"),
                                                          c("CMS1", "CMS4"),
                                                           c("CMS2", "CMS3"),
                                                           c("CMS2", "CMS4"),
                                                           c("CMS3", "CMS4")
                                                          )) +
                    #label = "p.signif") +
  ylab(gsub("_", " ", variable)) +
  xlab("") +
  ylim(0, 0.61)

violin_plot = plot + geom_violin(aes(fill = df_plot$CMS)) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = df_plot$CMS)) +
  ylab("") +
  theme(axis.text.x = element_blank())

pdf(file = "./Figures/TCR/013_Boxplot_clonality_CMS/Fig2c_Violin_plot_clonality_CMS_no_stats.pdf",
    width = 3, height = 3.8)
plot(violin_plot)
dev.off()

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

png(filename = "./Figures/TCR/013_Boxplot_clonality_CMS/Violin_plot_clonality_by_ICR_facet_by_CMS_stats.png",
    width = 6, height = 3, units = "in", res = 600)
plot(plot)
dev.off()
