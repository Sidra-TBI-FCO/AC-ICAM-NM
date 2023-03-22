
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
rownames(Genus_full_abundance)[grep("Ruminoc", rownames(Genus_full_abundance))]
Genus = "D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
  #"D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_01_Genus_by_ICR"), showWarnings = FALSE)

# Analysis
table_cluster_assignment = table_cluster_assignment[which(
  substring(rownames(table_cluster_assignment), 1, 3) %in% substring(colnames(Genus_full_abundance), 1, 3)),]

df_plot = data.frame(Patient_ID = substring(rownames(table_cluster_assignment), 1, 3),
                     ICR_cluster = table_cluster_assignment$ICR_HML,
                     ICRscore = table_cluster_assignment$ICRscore,
                     Abundance = NA)
df_plot$Abundance = Genus_full_abundance[Genus,][match(df_plot$Patient_ID,
                                                       substring(colnames(Genus_full_abundance), 1, 3))]

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$ICR_cluster = factor(df_plot$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

plot = ggplot(df_plot, aes(x = ICR_cluster, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, aes(fill=ICR_cluster)) +
  scale_fill_manual(values = c("ICR Low" = "blue", "ICR Medium" = "green", "ICR High" = "red")) +
  geom_jitter(width = 0.2, size = 0.8) +
  #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
  theme_bw() +
  scale_y_log10() +
  ylab("") +
  #ylab(paste0(gsub(".*\\D5__", "", Genus), "\n","Relative abundance")) +
  xlab("") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none")

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_01_Genus_by_ICR/EDF8a_", Genus, 
           "_by_ICR_cluster.pdf"), width = 3, height = 4)
plot(plot)
dev.off()

# run stats without running line 35  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
df_plot$ICR_cluster = as.numeric(df_plot$ICR_cluster)
cor.test(df_plot$Abundance, df_plot$ICR_cluster, method = "spearman")
