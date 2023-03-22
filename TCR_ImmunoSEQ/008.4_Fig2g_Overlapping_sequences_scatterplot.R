

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "gridExtra", "png", "grid")
ipak(required.packages)

Files = list.files("./Analysis/TCR/7_Overlapping_sequences_analysis/8_Per_patient_TN")
Patients = substring(Files, 1, 3)

Patients = Patients[2:10]

i = 9
for (i in 1:9){
  Patient = Patients[i]
  load(paste0("./Analysis/TCR/7_Overlapping_sequences_analysis/8_Per_patient_TN/", Patient, "_TN_clones_as_beausang.Rdata"))
  
  plot2 = ggplot(subset, aes(x = subset[,paste0(Patient, "N")], y = subset[,paste0(Patient, "T")])) +
    geom_point(aes(color = cat_beausang)) +
    scale_color_manual(values = c("red", "black", "black")) +
    scale_y_log10(limits = c(0.001, 10),breaks=c(0.01,0.1,1,10)) +
    scale_x_log10(limits = c(0.001, 10),breaks=c(0.01, 0.1,1,10)) +
    theme_bw() +
    ylab("") +
    xlab("") +
    #stat_bin2d(bins = 200) +
    #scale_fill_gradient(low = "black", high = "black", limits = c(0, 1000000000)) +
    coord_cartesian(clip = 'off') +
    theme(#axis.text.x = element_text(color = "black", size = 30),
          axis.text = element_blank(),
          #axis.text.y = element_text(color = "black", size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_blank(),
          plot.title = element_blank(),
          aspect.ratio = 1/1,
          legend.position = "none") +
    geom_segment(aes(x = 0.003, y = 0.10, xend = 0.5, yend = 10), linetype = "dashed") +
    geom_segment(aes(x = 0, y = 0.10, xend = 0.003, yend = 0.10), linetype = "dashed")
  
  png(paste0("./Figures/TCR/8.4_Enriched_sequences_scatterplots/0000_EDF4c_Without_label_Patient_", Patient, "_T_N_Scatterplot.png"), width = 7.3, height = 6,
      res=600, units = "in")
  plot(plot2)
  dev.off()
}
