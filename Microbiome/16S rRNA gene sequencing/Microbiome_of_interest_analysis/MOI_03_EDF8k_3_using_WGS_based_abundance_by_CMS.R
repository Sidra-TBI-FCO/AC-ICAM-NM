
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters 
dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_03_Genus_by_CMS"), showWarnings = FALSE)

# Set parameters 
rownames(Species_WGS)
Species = "s__Fusobacterium_nucleatum"
MSI = "" # or "MSI-H"

# Getting MSS and MSI-H
include_patients = MANTIS$Patient_ID[which(MANTIS$MSI == MSI)]

# Analysis
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"

df_plot = data.frame(Patient_ID = substring(colnames(Species_WGS), 1, 3),
                     CMS = NA,
                     Abundance = Species_WGS[Species,])

df_plot$CMS = Rfcms$RF.predictedCMS[match(df_plot$Patient_ID, substring(rownames(Rfcms), 1, 3))]

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$CMS = factor(df_plot$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
table(df_plot$CMS)

if(MSI == ""){}else{
  df_plot = df_plot[which(df_plot$Patient_ID %in% include_patients),]
}

plot = ggplot(df_plot, aes(x = CMS, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, aes(fill=CMS)) +
  scale_fill_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                               "mixed" = "lightgrey")) +
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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_03_Genus_by_CMS/EDF8k_3_WGS_based_", Species, "_", MSI,
           "_by_CMS.pdf"), width = 4, height = 4)
plot(plot)
dev.off()
 
# run stats without running line 34  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
t.test(df_plot$Abundance[which(df_plot$CMS == "CMS1")],
       df_plot$Abundance[-which(df_plot$CMS == "CMS1")])
