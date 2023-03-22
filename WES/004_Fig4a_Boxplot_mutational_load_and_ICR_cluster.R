
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "ggrepel")
ipak(required.packages)

# Set parameters

## Explore plot
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(file = paste0("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata"))
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]

Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"

df_plot = data.frame(Patient_ID = frequency_df$Patient_ID, Mutation_rate = frequency_df$Nonsilent_mutational_burden_per_Mb,
                     ICRscore = NA, HML_cluster = NA)
df_plot$Patient_ID = as.character(df_plot$Patient_ID)
df_plot$ICRscore = table_cluster_assignment$ICRscore[match(df_plot$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
df_plot$HML_cluster = table_cluster_assignment$ICR_HML[match(df_plot$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
df_plot$ESTIMATE_Immune = ESTIMATE[, c("ImmuneScore")][match(df_plot$Patient_ID, substring(rownames(ESTIMATE), 1, 3))]

df_plot$MANTIS = MANTIS$MSI[match(df_plot$Patient_ID, MANTIS$Patient_ID)]
df_plot$CMS = Rfcms$RF.predictedCMS[match(df_plot$Patient_ID, substring(rownames(Rfcms), 1, 3))]
df_plot$Mutation_frequency = frequency_df$Non_silent_Mutation_frequency
df_plot$histology = clinical_data$Tumor_morphology[match(df_plot$Patient_ID,
                                                         clinical_data$Patient_ID)]
                                                   

dir.create("./Figures/WES/004_Mutational_load_boxplots", showWarnings = FALSE)

colors = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
           "CMS3" = "#E97AA8", "CMS4" = "#009E74",
           "mixed" = "lightgrey")

MANTIS_colors = c("MSI-H" = "purple", "MSS" = "black")

histology_colors = c("adenocarcinoma nno" = "grey", "adenocarcinoom in villeus adenoom" = "green", 
                     "adenocarcinoom, intestinaal type" = "purple", "cribriform carcinoom" = "red", 
                     "mucineus adenocarcinoom" = "darkorange",
                     "zegelringcel carcinoom" = "darkgreen", "adenocarcinoom met gemengde subtypes" = "darkblue",
                     "neoplasma, maligne" = "pink")

plot = ggplot(df_plot, aes(y = Mutation_rate, x = ICRscore)) +
  #stat_cor(method = "pearson", size = 5) +
  geom_point(aes(color = CMS)) +
  scale_color_manual(values = colors) +
  scale_y_log10() +
  theme_bw() +
  #geom_text_repel() +
  ylab("Nonsynonymous Mutation load per Mb") +
  xlab("ICR score") +
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        aspect.ratio = 1/1)

png("./Figures/WES/004_Mutational_load_boxplots/Mutation_rate_by_ICRscore_scatter_CMS.png",
    width = 5.5, height = 5.5, units = "in", res = 600)
plot(plot)
dev.off()


df_plot$HML_cluster = factor(df_plot$HML_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

#84D976 = lightgreen (Mut. Low- ICR Low) 4
#7C7AED = lightpurple (Mut. High- ICR Low) 3
#00871F = darkgreen (Mut. Low- ICR High) 2
#8B16DD = darkpurple (Mut. High- ICR High) 1

# Remove ICR Medium
df_plot = df_plot[-which(df_plot$HML_cluster == "ICR Medium"),]
df_plot$Mutation_cat = NA
df_plot$Mutation_cat[which(df_plot$Mutation_rate < 12)] = "Mut_Low"
df_plot$Mutation_cat[which(df_plot$Mutation_rate >= 12)] = "Mut_High"
df_plot$ICR_Mut_Comb = paste(df_plot$HML_cluster, df_plot$Mutation_cat, sep = "")
table(df_plot$ICR_Mut_Comb)
plot = ggplot(df_plot, aes(y = Mutation_rate, x = HML_cluster)) +
  #scale_fill_manual(values = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.8, aes(color = ICR_Mut_Comb)) +
  scale_color_manual(values =  c("ICR LowMut_Low" = "#84D976", "ICR LowMut_High" = "#7C7AED",
                                 "ICR HighMut_Low" = "#00871F", "ICR HighMut_High" = "#8B16DD")) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE), limits = c(0.1, 400)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.position = "none") +
  #stat_compare_means(method = "t.test", comparisons = list(c("ICR High", "ICR Low")),
   #                  label = "p.signif") +
  ylab("Nonsynonymous mutation frequency \n (per Mb)") +
  xlab("")

pdf("./Figures/WES/004_Mutational_load_boxplots/Fig4a_Boxplot_Mutation_rate_by_ICR_cluster.pdf",
    width = 3.5, height = 4)
plot(plot)
dev.off()

test = t.test(df_plot$Mutation_frequency[which(df_plot$HML_cluster == "ICR Low")],
              df_plot$Mutation_frequency[which(df_plot$HML_cluster == "ICR High")])
test$p.value
