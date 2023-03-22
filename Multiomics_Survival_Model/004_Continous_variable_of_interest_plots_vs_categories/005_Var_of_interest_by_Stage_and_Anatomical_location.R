


# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Set parameters
var_of_interest = "var15_Microbiome_risk_score"
subset = "246"

# Load data
load("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")

if(subset == "246"){
  df = df[which(df$Microbiome_Paired_TN == "yes"),]
}

# Set parameters
dir.create(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots"), showWarnings = FALSE)

# Analysis
df_plot = df

df_plot$Anatomical_location = clinical_data$tumour_anatomic_site[match(df_plot$Patient_ID, clinical_data$Patient_ID)]

df_plot$Anatomical_location = factor(df_plot$Anatomical_location, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                             "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                             "rectosigmoideum"))

table(df_plot$Anatomical_location)

df_plot$var17_Stage_ordinal = factor(df_plot$var17_Stage_ordinal, levels = c("1", "2", "3", "4"))
table(df_plot$var17_Stage_ordinal)

plot = ggplot(df_plot, aes(x = var17_Stage_ordinal, y = get(var_of_interest))) +
  geom_boxplot(outlier.shape = NA, aes(fill = var17_Stage_ordinal)) +
  scale_fill_manual(values = c("1" = "#BFFFA4", "2" = "#FCFF8D", "3" = "#FFD580", "4" = "#FF6C6C")) +
  geom_jitter(width = 0.2, size = 0.8) +
  #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
  theme_bw() +
  ylab("") +
  #ylab(paste0(gsub(".*\\D5__", "", Genus), "\n","Relative abundance")) +
  xlab("") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none")

png(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots/", var_of_interest, "_by_Stage", 
           "_in_ICAM", subset, ".png"), res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()

plot = ggplot(df_plot, aes(x = Anatomical_location, y = get(var_of_interest))) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_manual(values = c("1" = "#BFFFA4", "2" = "#FCFF8D", "3" = "#FFD580", "4" = "#FF6C6C")) +
  geom_jitter(width = 0.2, size = 0.8) +
  #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
  theme_bw() +
  ylab("") +
  #ylab(paste0(gsub(".*\\D5__", "", Genus), "\n","Relative abundance")) +
  xlab("") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none")

png(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots/", var_of_interest, "_by_anatomical_location", 
           "_in_ICAM", subset, ".png"), res = 600, width = 8, height = 3.5, units = "in")
plot(plot)
dev.off()
