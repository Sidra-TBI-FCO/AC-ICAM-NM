

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters 
rownames(Species_WGS)
Species = "s__Fusobacterium_nucleatum"
MSI = "MSI-H" # or "MSI-H"

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_01_Genus_by_ICR"), showWarnings = FALSE)

# Getting MSS and MSI-H
include_patients = MANTIS$Patient_ID[which(MANTIS$MSI == MSI)]

# Analysis
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

table_cluster_assignment = table_cluster_assignment[which(
  substring(rownames(table_cluster_assignment), 1, 3) %in% substring(colnames(Species_WGS), 1, 3)),]

df_plot = data.frame(Patient_ID = substring(rownames(table_cluster_assignment), 1, 3),
                     ICR_cluster = table_cluster_assignment$ICR_HML,
                     ICRscore = table_cluster_assignment$ICRscore,
                     Abundance = NA)
df_plot$Abundance = Species_WGS[Species,][match(df_plot$Patient_ID,
                                                       substring(colnames(Species_WGS), 1, 3))]

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$ICR_cluster = factor(df_plot$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

if(MSI == ""){}else{
  df_plot = df_plot[which(df_plot$Patient_ID %in% include_patients),]
}

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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_01_Genus_by_ICR/EDF8h_WGS_based_", Species, "_", MSI,
           "_by_ICR_cluster.pdf"), width = 3, height = 4)
plot(plot)
dev.off()

# run stats without running line 35  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
df_plot$ICR_cluster = as.numeric(df_plot$ICR_cluster)
cor.test(df_plot$Abundance, df_plot$ICR_cluster, method = "spearman")
