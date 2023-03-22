

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Load data
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
MiXCR = read.csv("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Diversity.csv", 
                 stringsAsFactors = FALSE)
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
Exclude = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)

# Set parameters
Type = "TRB" # ALL IGH IGK IGL TRA TRB TRD TRG 

# Filter type
MiXCR$Type = substring(MiXCR$X, 8, 10)
MiXCR = MiXCR[which(MiXCR$Type == Type),]

# Prepare data
colnames(MiXCR)[1] = "Sample_ID"
MiXCR$Patient_ID = substring(MiXCR$Sample_ID, 1, 3)
MiXCR$Tissue = substring(MiXCR$Sample_ID, 4, 6)

# Filter only tumor samples
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),]

MiXCR$MiXCR_Clonality = 1- MiXCR$Pielou
MiXCR$ICRscore = table_cluster_assignment$ICRscore[match(substring(MiXCR$Sample_ID, 1, 4), 
                                                         substring(rownames(table_cluster_assignment), 1, 4))]
MiXCR$ICR_cluster = table_cluster_assignment$ICR_HML[match(substring(MiXCR$Sample_ID, 1, 4), 
                                                         substring(rownames(table_cluster_assignment), 1, 4))]

Exclude$Patient_ID = str_pad(Exclude$Patient_ID, pad = "0", 3)
MiXCR = MiXCR[-which(MiXCR$Patient_ID %in% Exclude$Patient_ID),]

save(MiXCR, file = "./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Clean_TRB_clonality_341_patients.Rdata")

#MiXCR = MiXCR[which(MiXCR$Cohort == "validation"),]

df_plot = MiXCR[, c("MiXCR_Clonality", "ICRscore", "ICR_cluster", "Sample_ID")]
df_plot$ICR_cluster = factor(df_plot$ICR_cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))

plot = ggplot(df_plot, aes(x = ICRscore, y = MiXCR_Clonality)) +
  geom_point(aes(color = ICR_cluster), size = 0.8) +
  ylab("MiXCR Clonality") +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        aspect.ratio = 1/1,
        legend.position = "none") +
  labs(color = "ICR cluster") +
  #stat_cor(method = "pearson", size = 5, digits = 2) +
  geom_smooth(method="lm")

dir.create("./Figures/TCR/019_MiXCR_ICRscore_in_validation_cohort", showWarnings = FALSE)
pdf(paste0("./Figures/TCR/019_MiXCR_ICRscore_in_validation_cohort/019_Fig2b_Scatter_MiXCR_ICRscore_Complete_cohort.pdf"),
    width = 3, height = 3)
plot(plot)
dev.off()
