
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

# Load data 

# Set parameters
clones = "Tumor restricted"         # which clones (Tumor restricted, Overlapping, Normal restricted)
tissue = "T"                   # frequency in which sample (T or N)
Files = list.files("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN")
Patients = substring(Files, 1, 3)

dir.create("./Figures/TCR/7_Overlapping_sequences_scatterplots", showWarnings = FALSE)

if(clones == "Normal restricted"){
  color = "#F57900"
}
if(clones == "Overlapping"){
  color = "#00BA38"
}
if(clones == "Tumor restricted"){
  color = "#619CFF"
}

if(tissue == "N"){
  text = " T-cell clones \n in normal colon tissue"
}
if(tissue == "T"){
  text = " T-cell clones \n in tumor tissue"
}

i=2
for (i in 1:10){
  Patient = Patients[i]
  load(paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/7.2_Per_patient_TN/", Patient, "_TN_clones.Rdata"))
  
  subset = subset[which(subset$cat == clones),]
  
  plot = qplot(subset[,paste0(Patient, tissue)],
               bins = 60,
               main = paste0("Patient ", Patient, text),
               xlab = "Productive Frequency (%)",
               ylab = "frequency (n of clones)",
               fill = I(color),
               col = I("black")) +
    scale_x_log10(limits = c(0.001, 10),breaks=c(0.01, 0.1,1,10)) +
    scale_y_log10() +
    theme_bw() +
    theme(axis.text.x = element_text(size=15, color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title = element_text(size=15, color = "black"),
          title = element_text(size=15, face="bold"),
          aspect.ratio = 1/3.5) 
    
  
  dir.create("./Figures/TCR/7_Overlapping_sequences_scatterplots/7.2_Histograms_Clone_distribution", showWarnings =FALSE)
  pdf(paste0("./Figures/TCR/7_Overlapping_sequences_scatterplots/7.2_Histograms_Clone_distribution/7.2_EDF4b_Patient_", Patient, "_Histogram_",  
             clones, "_", tissue,".pdf"), width = 5, height = 3)
  print(plot)
  dev.off()
  
  plot = qplot(subset[,paste0(Patient, tissue)],
               bins = 60,
               main = paste0("Patient ", Patient, text),
               xlab = "",
               ylab = "",
               fill = I(color),
               col = I("black")) +
    scale_x_log10(limits = c(0.001, 10),breaks=c(0.01, 0.1,1,10)) +
    scale_y_log10() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          title = element_blank(),
          aspect.ratio = 1/3.5) 
  dir.create("./Figures/TCR/7_Overlapping_sequences_scatterplots/7.2_Histograms_Clone_distribution_without_labels", showWarnings =FALSE)
  png(paste0("./Figures/TCR/7_Overlapping_sequences_scatterplots/7.2_Histograms_Clone_distribution_without_labels/7.2_Patient_", Patient, "_Histogram_",  
             clones, "_", tissue,".png"), width = 5, height = 3,
      units = "in", res = 600)
  print(plot)
  dev.off()
  
}

