
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Load data
Best_CV_coef = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv", stringsAsFactors = FALSE)

colnames(Best_CV_coef) = c("rownames", colnames(Best_CV_coef)[1:3])
Best_CV_coef$Family_Genus = gsub(".*\\D_4__", "", Best_CV_coef$covariates)
Best_CV_coef$Family_Genus= gsub("D_5__", " - ", Best_CV_coef$Family_Genus)

Best_CV_coef$coefficients = as.numeric(Best_CV_coef$coefficients)
Best_CV_coef = Best_CV_coef[order(Best_CV_coef$coefficients, decreasing = TRUE),]
Best_CV_coef$Family_Genus = factor(Best_CV_coef$Family_Genus, levels = rev(Best_CV_coef$Family_Genus))

plot = ggplot(Best_CV_coef, aes(x = coefficients, y = Family_Genus, fill = colors)) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("") +
  theme(axis.text.x = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14),
        legend.position = "none")

dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/103_Coefficients_Risk_scores", showWarnings = FALSE)
png(file = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/103_Coefficients_Risk_scores/version2_Coefficients_Risk_scores.png"), 
    res = 600, units = "in", width = 8, height = 10)
plot(plot)
dev.off()

pdf(file = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/103_Coefficients_Risk_scores/version2_Coefficients_Risk_scores.pdf"), 
    width = 8, height = 10)
plot(plot)
dev.off()

