
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

# Set parameters
variable = "productive_frequency_sum" # "productive_frequency_sum" #"Fraction_Tumor_Enriched"

# Load data
load("./Analysis/TCR/8_As_Beausang_TN_Enriched/8_Tumor_enriched_sequences.Rdata")
if(variable == "productive_frequency_sum"){
  load("./Analysis/TCR/8_As_Beausang_TN_Enriched/8.3_Tumor_enriched_sequences_stats.Rdata")
}
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")

plot_df = data.frame(Patient = results$Patient, Variable = results[, variable], ICRscore = NA, ICR_cluster = NA)
plot_df$ICRscore = table_cluster_assignment$ICRscore[match(plot_df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
plot_df$ICR_cluster = table_cluster_assignment$ICR_HML[match(plot_df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
plot_df = plot_df[-which(plot_df$Patient == "175"),]

plot = ggplot(plot_df, aes(x = plot_df$ICRscore, y = plot_df$Variable, label = plot_df$Patient)) +
  geom_point(aes(color = plot_df$ICR_cluster)) +
  scale_color_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  #geom_abline() +
  ylab(gsub("_", " ", variable)) +
  xlab("ICR score") +
  theme_bw() +
  stat_cor(method = "pearson", size = 10) +
  geom_smooth(method="lm", se = FALSE, color = "black", size = 0.6) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 23, color = "black"),
        axis.text.y = element_text(size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"),
        aspect.ratio = 1/1) +
  geom_text(hjust = 0, nudge_x = 0.05, size = 7) +
  xlim(3.5, 8.5) 
  #xlab("") +
  #ylab("")

if(variable == "productive_frequency_sum"){
  plot = plot + ylab("Productive frequency \ntumor enriched total (%)")
  dir.create("./Figures/TCR/8.3_Scatterplot_Tumor_enriched_fraction", showWarnings = FALSE)
  png(paste0("./Figures/TCR/8.3_Scatterplot_Tumor_enriched_fraction/Dec_2022_Scatterplot_", variable, "_as_Beausang_by_ICRscore.png"),
      width = 5, height = 5, units = "in", res = 600)
  plot(plot)
  dev.off()
  
  pdf(paste0("./Figures/TCR/8.3_Scatterplot_Tumor_enriched_fraction/Dec_2022_Scatterplot_", variable, "_as_Beausang_by_ICRscore.pdf"),
      width = 5, height = 5)
  plot(plot)
  dev.off()
}else{
  dir.create("./Figures/TCR/8.1_Scatterplot_Tumor_enriched_fraction", showWarnings = FALSE)
  png(paste0("./Figures/TCR/8.1_Scatterplot_Tumor_enriched_fraction/v2_Scatterplot_", variable, "_as_Beausang_by_ICRscore.png"),
      width = 5, height = 5, units = "in", res = 600)
  plot(plot)
  dev.off()
}



