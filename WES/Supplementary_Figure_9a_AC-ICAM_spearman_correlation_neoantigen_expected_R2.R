# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db

setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))


required.bioconductor.packages = c("ggplot2", "ggpubr", "dplyr")                                                                   
ibiopak(required.bioconductor.packages)

#load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")

#colnames(Merged_dataset)

plot = ggplot(Merged_dataset, aes(x= Expected_neoantigen_count, y = Neoantigen_count, fill = Mutation_cat)) + 
  geom_point() +
  scale_color_manual(values = c("#EAAED0", "#A7EABD")) +
  stat_cor(method = "spearman", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Expected Neoantigen Count") +
  ylab("Neoantigen Count") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

print(plot)

# See script 020_Neoantigen_load_calculation for analogous AC-ICAM script
png("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_R2_colored.png", res = 600, width = 5.2, height = 5, units = "in") 
plot(plot)
dev.off()


png(paste0("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_colored.png"), res = 600, width = 4, height = 4, units = "in")
ggscatter(
  Merged_dataset,
  x= "Expected_neoantigen_count", y= "Neoantigen_count",
  fill = "Mutation_cat",
  shape = 21,
  size = 2,
  add = "reg.line",
  cor.coef = TRUE,
  cor.method = "spearman",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE,
  title=" ",
  xlab="Expected Neoantigen Count",
  ylab = "Neoantigen Count",
  palette = c("#EAAED0", "#A7EABD"))
dev.off()

#plot(plot)
#dev.off()

# supplementary figure 9a
png(paste0("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_slope_1_intercept_0.png"), 
    res = 600, width = 4, height = 4, units = "in")
plot = ggplot(Merged_dataset, aes(x= Expected_neoantigen_count, y = Neoantigen_count)) + 
  geom_point() +
  #stat_cor(method = "spearman", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Expected Neoantigen Count") +
  ylab("Observed Neoantigen Count") +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        aspect.ratio = 1/1)
  #geom_abline(intercept = 0, slope =1) +
  #geom_smooth(method="lm")
print(plot)
dev.off()

pdf(paste0("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_slope_1_intercept_0.pdf"), 
    width = 4, height = 4)
print(plot)
dev.off()



png(paste0("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_slope_1.png"), res = 600, width = 4, height = 4, units = "in")
plot = ggplot(Merged_dataset, aes(x= Expected_neoantigen_count, y = Neoantigen_count)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Expected Neoantigen Count") +
  ylab("Neoantigen Count") +
  theme(axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        aspect.ratio = 1/1) +
  geom_abline(slope =1) +
  geom_smooth(method="lm")
print(plot)
dev.off()

## remove regression line
png(paste0("./Microbiome_EA/Figures/spearman_correlation_observed_expected_neoantigen_slope_1_no_regression.png"), res = 600, width = 4, height = 4, units = "in")
plot = ggplot(Merged_dataset, aes(x= Expected_neoantigen_count, y = Neoantigen_count)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Expected Neoantigen Count") +
  ylab("Neoantigen Count") +
  theme(axis.title.x = element_text(size = 8, colour = "black"),
        axis.title.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        aspect.ratio = 1/1) +
  geom_abline(intercept = 0, slope =1)
  #geom_smooth(method="lm")
print(plot)
dev.off()

################################################################################
################################################################################

TMB = Merged_dataset

mutation = "hypermutated" # nonhypermutated # hypermutated

TMB = TMB[which(TMB$Mutation_cat == mutation),]

png(paste0("./Microbiome_EA/Figures/spearman_correlation_TMB_GIE_",mutation,".png"), res = 600, width = 4, height = 4, units = "in")
plot = ggplot(TMB, aes(x= Immunoediting_score, y = Nonsilent_mutational_burden_per_Mb)) + 
  geom_point() +
  stat_cor(method = "spearman", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Immunoediting_score") +
  ylab("Nonsilent_mutational_burden_per_Mb") +
  theme(axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

print(plot)
dev.off()

TMB$Immunoediting_score

################################################################################
################################################################################

# Boxplot

png(paste0("./Microbiome_EA/Figures/boxplot_TMB_GIE_categories.png"), res = 600, width = 2, height = 4, units = "in")

boxplot_gene = ggplot(data = Merged_dataset, aes(x= Immunoedited, y= Nonsilent_mutational_burden_per_Mb, fill = Immunoedited)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = c("#FF8C00", "#000080")) +
  scale_y_continuous("Nonsilent_mutational_burden_per_Mb") +
  ggtitle(" ") +
  stat_compare_means(method = "t.test") +
  geom_point(position=position_jitterdodge(),alpha=0.3, size = 0.05) + 
  theme(plot.title = element_text(size=8, face = "bold")) + theme(axis.title = element_text(size = 8, face = "bold", colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(angle = 90, colour = "black", size = 8)) +
  theme(axis.line = element_line(color= "black", size = 0.4)) +
  guides(fill=FALSE)

print(boxplot_gene)
dev.off()


















