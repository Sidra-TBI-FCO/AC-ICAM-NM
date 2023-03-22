

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

dir.create("./Figures/Exploration_reviewers/ES_by_anatomical_location", showWarnings = FALSE)
dir.create("./Analysis/Exploration_reviewers/ES_by_anatomical_location", showWarnings = FALSE)

ipak(c("RColorBrewer", "ggplot2", "easyGgplot2", "stringr"))

# Set parameters
Risk_score = "Tumor" # "Tumor" or "Normal"

# Load data
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

if(Risk_score == "Tumor"){
  load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
}

if(Risk_score == "Normal"){
  Training = read.csv("./Processed_Data/Microbiome/External/9_August/Normal_Risk_Scores.csv")
  Training$Patient_ID = str_pad(Training$samples, pad = "0", 3)
}


# Spearman correlation anatomical location and enrichment score
df = data.frame(Patient_ID = clinical_data$Patient_ID, Anatomic_location = clinical_data$tumour_anatomic_site,
                MBR = NA)
df$MBR = Training$prediction[match(df$Patient_ID, Training$Patient_ID)]
df = df[-which(is.na(df$MBR)),]

df$Anatomic_location[which(df$Anatomic_location == "ceceum")] = "caecum"

df$Anatomic_location = factor(df$Anatomic_location, levels = c("caecum", "colon ascendens", "flexura hepatica", 
                                                               "colon transversum", "flexura lienalis", "colon descendens", 
                                                               "colon sigmoideum", "rectosigmoideum"))
table(df$Anatomic_location)
df$Anatomic_location = as.numeric(df$Anatomic_location)

cor = cor.test(df$MBR, df$Anatomic_location, method = "spearman")
cor$estimate
cor$p.value

df_plot = df  
df_plot$Anatomic_location = factor(df_plot$Anatomic_location)
levels(df_plot$Anatomic_location) = c("caecum", "colon ascendens", "flexura hepatica", 
                                      "colon transversum", "flexura lienalis", "colon descendens", 
                                      "colon sigmoideum", "rectosigmoideum")
  
plot = ggplot(df_plot, aes(x = Anatomic_location, y = MBR, fill = Anatomic_location)) +
  scale_fill_manual(values = brewer.pal(n = 8, name = "Set3")) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.8) +
  ylab("MBR") +
  theme_bw() +
  xlab("") +
  theme(axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black",
                                   angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none") +
  geom_text(x = 5, y = max(df_plot$MBR), label = paste0("Rho = ", round(cor$estimate, 2), ", ", "p  = ", 
                                                       formatC(cor$p.value, format = "e", digits = 2)),
            check_overlap = TRUE)
  
pdf(paste0("./Figures/Exploration_reviewers/ES_by_anatomical_location/Boxplot_MBR_", Risk_score, "_by_anatomic_location.pdf"),
    width = 5, height = 4)
plot(plot)
dev.off()

