
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Load data
load("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")

# Prepare data for plot
df$var17_Stage_ordinal = factor(df$var17_Stage_ordinal)

plot = ggplot(df, aes(x = var17_Stage_ordinal, y = var15_Microbiome_risk_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.2) +
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"),
                                        c("3", "4")), 
                     method = "t.test",
                     label = "p.signif") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14))

dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/102_Risk_scores_by_stage", showWarnings = FALSE)

png("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/102_Risk_scores_by_stage/102_Risk_scores_by_stage.png",
    res = 600, units = "in", width = 6, height = 4)
plot(plot)
dev.off()

df$var02_ICR_cluster = factor(df$var02_ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
plot = ggplot(df, aes(x = var02_ICR_cluster, y = var15_Microbiome_risk_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.2) +
  stat_compare_means(comparisons = list(c("ICR Low", "ICR Medium"), c("ICR Medium", "ICR High"),
                                        c("ICR Low", "ICR High")), 
                     method = "t.test",
                     label = "p.signif") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14))

png("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/102_Risk_scores_by_stage/102_Risk_scores_by_ICR_cluster.png",
    res = 600, units = "in", width = 6, height = 4)
plot(plot)
dev.off()

df$var12_Hypermutation_status
plot = ggplot(df, aes(x = var12_Hypermutation_status, y = var15_Microbiome_risk_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.2) +
  stat_compare_means(comparisons = list(c("nonhypermutated", "hypermutated")), 
                     method = "t.test",
                     label = "p.signif") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14))


png("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/102_Risk_scores_by_stage/102_Risk_scores_by_Hypermutation_status.png",
    res = 600, units = "in", width = 6, height = 4)
plot(plot)
dev.off()

plot = ggplot(df, aes(x = var19_MSI_MANTIS, y = var15_Microbiome_risk_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, width = 0.2) +
  stat_compare_means(comparisons = list(c("MSS", "MSI-H")), 
                     method = "t.test",
                     label = "p.signif") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14))

png("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/102_Risk_scores_by_stage/102_Risk_scores_by_MSI_status.png",
    res = 600, units = "in", width = 6, height = 4)
plot(plot)
dev.off()
