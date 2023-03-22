
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

# Load data 
CombinedRearrangements = read.csv("./Processed_Data/TCR/CombinedRearrangements_v1.tsv", sep = "\t")
colnames(CombinedRearrangements) = gsub("X", "", colnames(CombinedRearrangements))

Normal = colnames(CombinedRearrangements)[grep("N", colnames(CombinedRearrangements))]
Matched_T = gsub("N", "T", Normal)
Patients = gsub("N", "", Normal)

# TN pairs
TN_Rearrangements = CombinedRearrangements[, c("Amino.Acid", "Sum..Productive.Frequency.", "Present.In", Normal, Matched_T)]

dir.create("./Figures/TCR/7_Overlapping_sequences_scatterplots", showWarnings = FALSE)

results = data.frame(Patient = Patients, Number_Tumor_Restricted = NA, Number_Normal_Restricted = NA, Number_Overlap = NA,
                     Median_PF_Tumor = NA, Median_PF_Normal = NA, third_quartile_T = NA, third_quartile_N = NA)

i=2
for (i in 1:9){
  Patient = Patients[i]
  subset = TN_Rearrangements[, c("Amino.Acid", "Sum..Productive.Frequency.", paste0(Patient, "N"), paste0(Patient, "T"))]
  subset = subset[-which(subset[,paste0(Patient, "N")] == 0 & subset[,paste0(Patient, "T")] == 0),]
  subset[,paste0(Patient, "N")] = subset[,paste0(Patient, "N")] *100
  subset[,paste0(Patient, "T")] = subset[,paste0(Patient, "T")] *100
  subset$cat = "Overlapping"
  subset$cat[which(subset[,paste0(Patient, "N")] == 0)] = "Tumor restricted"
  subset$cat[which(subset[,paste0(Patient, "T")] == 0)] = "Normal restricted"
  plot = ggplot(subset, aes(x = subset[,paste0(Patient, "N")], y = subset[,paste0(Patient, "T")])) +
    geom_point(aes(color = cat)) +
    scale_color_manual(labels = c(paste0("Normal restricted (n=", length(which(subset$cat == "Normal restricted")),")"),
                                  paste0("Overlapping (n=", length(which(subset$cat == "Overlapping")), ")"),
                                  paste0("Tumor restricted (n=", length(which(subset$cat == "Tumor restricted")), ")")),
                       values = c("#F57900", "#00BA38", "#619CFF")) +
    scale_y_log10(limits = c(0.001, 10),breaks=c(0.01,0.1,1,10)) +
    scale_x_log10(limits = c(0.001, 10),breaks=c(0.01, 0.1,1,10)) +
    theme_bw() +
    ylab(paste0("Tumor SUM (Productive Frequency)")) +
    xlab(paste0("Normal SUM (Productive Frequency)")) +
    coord_cartesian(clip = 'off') +
    theme(axis.text.x = element_text(color = "black", size = 20),
          axis.text.y = element_text(color = "black", size = 20),
          axis.title.x = element_text(color = "black", size = 20),
          axis.title.y = element_text(color = "black", size = 20),
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size = 14),
          plot.title = element_text(color = "black", face = "bold", size = 22),
          aspect.ratio = 1/1,
          legend.position = "bottom") +
    ggtitle(paste0("Patient ", Patient))
  
  plot2 = ggplot(subset, aes(x = subset[,paste0(Patient, "N")], y = subset[,paste0(Patient, "T")])) +
    geom_point(aes(color = cat)) +
    scale_color_manual(values = c("#F57900", "#00BA38", "#619CFF")) +
    scale_y_log10(limits = c(0.001, 10),breaks=c(0.01,0.1,1,10)) +
    scale_x_log10(limits = c(0.001, 10),breaks=c(0.01, 0.1,1,10)) +
    theme_bw() +
    ylab("") +
    xlab("") +
    coord_cartesian(clip = 'off') +
    theme(axis.text.x = element_text(color = "black", size = 20),
          axis.text.y = element_text(color = "black", size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_blank(),
          plot.title = element_blank(),
          aspect.ratio = 1/1,
          legend.position = "none")
  
  pdf(paste0("./Figures/TCR/7_Overlapping_sequences_scatterplots/7_EDF4b_Patient_", Patient, "_T_N_Scatterplot.pdf"), width = 9, height = 6)
  plot(plot2)
  dev.off()
  
  dir.create(paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN"), showWarnings = FALSE)
  save(subset, file = paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN/", Patient, "_TN_clones.Rdata"))
  
  table = table(subset$cat)
  results$Number_Tumor_Restricted[which(results$Patient == Patient)] = table["Tumor restricted"]
  results$Number_Normal_Restricted[which(results$Patient == Patient)] = table["Normal restricted"]
  results$Number_Overlap[which(results$Patient == Patient)] = table["Overlapping"]
  PF_clones_T = subset[,paste0(Patient, "T")][-which(subset[,paste0(Patient, "T")] == 0)]
  PF_clones_N = subset[,paste0(Patient, "N")][-which(subset[,paste0(Patient, "N")] == 0)]
  results$Median_PF_Tumor[which(results$Patient == Patient)] = median(PF_clones_T)
  results$Median_PF_Normal[which(results$Patient == Patient)] = median(PF_clones_N)
  results$third_quartile_T[which(results$Patient == Patient)] = quantile(PF_clones_T, 3/4)
  results$third_quartile_N[which(results$Patient == Patient)] = quantile(PF_clones_N, 3/4)
}
