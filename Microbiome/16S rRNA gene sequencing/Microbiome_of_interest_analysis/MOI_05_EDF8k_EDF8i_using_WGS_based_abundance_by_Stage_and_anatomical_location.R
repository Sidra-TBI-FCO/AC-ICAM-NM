
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters 
rownames(Species_WGS)
Species = "s__Fusobacterium_nucleatum"
MSI = "" # or "MSI-H"

# Getting MSS and MSI-H
include_patients = MANTIS$Patient_ID[which(MANTIS$MSI == MSI)]

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage"), showWarnings = FALSE)

# Analysis
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

df_plot = data.frame(Patient_ID = substring(colnames(Species_WGS), 1, 3),
                     Stage = NA,
                     Anatomical_location = NA,
                     Abundance = Species_WGS[Species,])

df_plot$Stage = clinical_data$ajcc_pathologic_tumor_stage[match(df_plot$Patient_ID, clinical_data$Patient_ID)]
df_plot$Anatomical_location = clinical_data$tumour_anatomic_site[match(df_plot$Patient_ID, clinical_data$Patient_ID)]

df_plot$Anatomical_location = factor(df_plot$Anatomical_location, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                             "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                             "rectosigmoideum"))

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$Stage = factor(df_plot$Stage, levels = c("1", "2", "3", "4"))
table(df_plot$Stage)

if(MSI == ""){}else{
  df_plot = df_plot[which(df_plot$Patient_ID %in% include_patients),]
}


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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage/EDF8k_4_WGS_based_", Species, "_", MSI,
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_05_Genus_by_Stage/EDF8i_WGS_based_", Species, "_", MSI, 
           "_by_Anatomical_location.pdf"), width = 8, height = 3.5)
plot(plot)
dev.off()

# run stats without running line 39  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
df_plot$Stage = as.numeric(df_plot$Stage)
cor.test(df_plot$Abundance, df_plot$Stage)

