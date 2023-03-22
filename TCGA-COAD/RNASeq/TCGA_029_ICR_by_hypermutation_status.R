
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("ggplot2", "ggpubr")                                                                   
ibiopak(required.bioconductor.packages)

# Load data
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Deconvolution_and_GSEA/ICR_genes_ES.Rdata")

# Set parameters
ICR = "ICRscore" # "ICRscore" or "ICR_ES"

# Analysis
plot_df = frequency_df

plot_df$Mutation_cat = NA
plot_df$Mutation_cat[which(plot_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
plot_df$Mutation_cat[which(plot_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

plot_df$Mutation_cat = factor(plot_df$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))

plot_df$ICRscore = table_cluster_assignment$ICRscore[match(plot_df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 12))]
plot_df$ICR_ES = ES[1,][match(plot_df$Patient_ID, substring(colnames(ES), 1, 12))]

plot_df = plot_df[-which(is.na(plot_df$ICR_ES)),]

plot = ggplot(plot_df, aes(x = Mutation_cat, y = get(ICR), fill = Mutation_cat)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("nonhypermutated" = "#A7EABD", "hypermutated" = "#EAAED0")) +
  geom_jitter(width = 0.2, size = 0.8) +
  theme_bw() +
  xlab("") +
  ylim(3.8, 11.5) +
  ylab(paste0(ICR)) +
  stat_compare_means(method = "t.test", comparisons = list(c("nonhypermutated", "hypermutated"))) +
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        title = element_text(size = 15, color = "black"),
        legend.position = "none")

dir.create("./Figures/029_ICR_by_hypermutation_status", showWarnings = FALSE)
png(paste0("./Figures/029_ICR_by_hypermutation_status/", ICR, "_by_hypermutation_status.png"),
    res = 600, units = "in", width = 3, height = 5.5)
plot(plot)
dev.off()
