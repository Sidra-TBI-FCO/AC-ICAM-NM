
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("reshape2", "ggplot2")
ipak(required.packages)

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Raw/JSREP.Complete.gene.filtered.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")

expression.filtered.log = log(expression.filtered +1, 2)

dir.create("./Figures/Trimmed_p", showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Normalization_Boxplots", showWarnings = FALSE)
png("./Figures/Trimmed_p/Normalization_Boxplots/Pre-Normalization_Boxplot.png", res = 600, height = 5, width = 10, units = "in")
boxplot(expression.filtered.log, las = 2, ylab = "Log2 transformed expression values")
dev.off()

png("./Figures/Trimmed_p/Normalization_Boxplots/Post-Normalization_Boxplot.png", res = 600, height = 5, width = 10, units = "in")
boxplot(RNASeq.QN.LOG2, las = 2, ylab = "Log2 transformed normalized expression values")
dev.off()

df_plot_pre = melt(expression.filtered.log)
colnames(df_plot_pre) = c("Gene", "Sample", "Count")
png("./Figures/Trimmed_p/Normalization_Boxplots/Pre-Normalization_Density_plot.png", res = 600, height = 5, width = 10, units = "in")
plot = ggplot(df_plot_pre, aes(x = Count)) +
  geom_density(aes(group = Sample), show.legend = NA) +
  theme_bw()
plot(plot)
dev.off()

df_plot_post = melt(RNASeq.QN.LOG2)
colnames(df_plot_post) = c("Gene", "Sample", "Count")
png("./Figures/Trimmed_p/Normalization_Boxplots/Post-Normalization_Density_plot.png", res = 600, height = 5, width = 10, units = "in")
plot = ggplot(df_plot_post, aes(x = Count)) +
  geom_density(aes(group = Sample)) +
  theme_bw()
plot(plot)
dev.off()
