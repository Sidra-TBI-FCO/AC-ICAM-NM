

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "stringr", "survminer", "survival", "coin", "survival", "dplyr"))

dir.create("./Figures/Multiomics_Survival_Model/004_Pie_Charts_Stage", showWarnings = FALSE)

# Load data
load("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
Merged_dataset = merge(df, clinical_data, by = "Patient_ID")

Merged_dataset$var02_ICR_cluster = factor(Merged_dataset$var02_ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$Combined_group = paste(Merged_dataset$var02_ICR_cluster, Merged_dataset$var16_Microbiome_risk_group, sep = "- ")
Merged_dataset$Combined_group[-which(Merged_dataset$Combined_group == "ICR High- Low Risk")] = "Rest"
Merged_dataset$Combined_group = factor(Merged_dataset$Combined_group, levels = c("ICR High- Low Risk", "Rest"))

# Chi-squared test
table(Merged_dataset$Combined_group, Merged_dataset$var17_Stage_ordinal)
chisq.test(table(Merged_dataset$Combined_group, Merged_dataset$var17_Stage_ordinal))

Merged_dataset = Merged_dataset[which(Merged_dataset$Combined_group == "ICR High- Low Risk"),]
Merged_dataset = Merged_dataset[which(Merged_dataset$Microbiome_Paired_TN == "yes"),]

DF1 <- Merged_dataset %>%
  group_by(ajcc_pathologic_tumor_stage) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

colors = c(c("1" = "#55B78D", "2" = "#F3EB82",
             "3" = "#F3C182", "4" = "#F66258"))

plot = ggplot(DF1, aes(x = "", y =perc*100, fill = ajcc_pathologic_tumor_stage)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "ajcc_pathologic_tumor_stage", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors) + coord_polar("y", start=0)


png(paste0("./Figures/Multiomics_Survival_Model/004_Pie_Charts_Stage/Stage_Pie_Chart_AC_ICAM_in_ICR_High_MBR_Low.png"),
    height = 7, width = 7, units = "in", res = 600)
plot(plot)
dev.off()

pdf(paste0("./Figures/Multiomics_Survival_Model/004_Pie_Charts_Stage/Stage_Pie_Chart_AC_ICAM_in_ICR_High_MBR_Low.pdf"),
    height = 7, width = 7)
plot(plot)
dev.off()






 