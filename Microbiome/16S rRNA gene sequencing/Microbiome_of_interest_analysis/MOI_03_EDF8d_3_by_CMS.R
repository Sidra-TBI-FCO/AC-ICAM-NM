
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
# rownames(Genus_full_abundance)[grep("Fuso", rownames(Genus_full_abundance))
Genus = "D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
  #"D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 1"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_03_Genus_by_CMS"), showWarnings = FALSE)

# Analysis
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"

df_plot = data.frame(Patient_ID = substring(colnames(Genus_full_abundance), 1, 3),
                     CMS = NA,
                     Abundance = Genus_full_abundance[Genus,])

df_plot$CMS = Rfcms$RF.predictedCMS[match(df_plot$Patient_ID, substring(rownames(Rfcms), 1, 3))]

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$CMS = factor(df_plot$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
table(df_plot$CMS)

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

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_03_Genus_by_CMS/EDF8d_3_", Genus, 
           "_by_CMS.pdf"), width = 4, height = 4)
plot(plot)
dev.off()

t.test(df_plot$Abundance[which(df_plot$CMS == "CMS1")],
       df_plot$Abundance[-which(df_plot$CMS == "CMS1")])
