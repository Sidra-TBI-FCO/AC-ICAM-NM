

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
# rownames(Genus_full_abundance)[grep("Fuso", rownames(Genus_full_abundance))
Genus = "D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
  # "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 1"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage"), showWarnings = FALSE)

# Analysis

df_plot = data.frame(Patient_ID = substring(colnames(Genus_full_abundance), 1, 3),
                     Stage = NA,
                     Anatomical_location = NA,
                     Abundance = Genus_full_abundance[Genus,])

df_plot$Stage = clinical_data$ajcc_pathologic_tumor_stage[match(df_plot$Patient_ID, clinical_data$Patient_ID)]
df_plot$Anatomical_location = clinical_data$tumour_anatomic_site[match(df_plot$Patient_ID, clinical_data$Patient_ID)]

df_plot$Anatomical_location = factor(df_plot$Anatomical_location, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                             "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                             "rectosigmoideum"))

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$Stage = factor(df_plot$Stage, levels = c("1", "2", "3", "4"))
table(df_plot$Stage)

plot = ggplot(df_plot, aes(x = Stage, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Stage)) +
  scale_fill_manual(values = c("1" = "#BFFFA4", "2" = "#FCFF8D", "3" = "#FFD580", "4" = "#FF6C6C")) +
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage/EDF8d_", Genus, 
           "_by_Stage.pdf"), width = 4, height = 4)
plot(plot)
dev.off()

plot = ggplot(df_plot, aes(x = Anatomical_location, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_manual(values = c("1" = "#BFFFA4", "2" = "#FCFF8D", "3" = "#FFD580", "4" = "#FF6C6C")) +
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage/EDF8b_", Genus, 
           "_by_Anatomical_location.pdf"), width = 8, height = 3.5)
plot(plot)
dev.off()

df_plot$Stage = as.numeric(df_plot$Stage)
cor.test(df_plot$Abundance, df_plot$Stage)

