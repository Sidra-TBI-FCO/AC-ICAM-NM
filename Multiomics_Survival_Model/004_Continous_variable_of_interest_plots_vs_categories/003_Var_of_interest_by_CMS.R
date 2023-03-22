


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

if(subset == "246"){
  df = df[which(df$Microbiome_Paired_TN == "yes"),]
}

# Set parameters
dir.create(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots"), showWarnings = FALSE)

# Analysis
df_plot = df
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
df_plot$CMS = Rfcms$RF.predictedCMS[match(df_plot$Patient_ID, substring(rownames(Rfcms), 1, 3))]

plot = ggplot(df_plot, aes(x = CMS, y = get(var_of_interest))) +
  geom_boxplot(outlier.shape = NA, aes(fill= CMS)) +
  scale_fill_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                               "mixed" = "lightgrey")) +
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

png(paste0("./Figures/Multiomics_Survival_Model/004_Var_of_interest_plots/", var_of_interest, 
           "_by_CMS_in_ACICAM", subset, ".png"), res = 600, width = 4, height = 4, units = "in")
plot(plot)
dev.off()
