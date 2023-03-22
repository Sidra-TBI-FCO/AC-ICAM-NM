

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "stringr", "survminer", "survival", "coin", "survival"))

dir.create("./Figures/Multiomics_Survival_Model", showWarnings = FALSE)
dir.create("./Figures/Multiomics_Survival_Model/002_Multivariate_survival_analysis", showWarnings = FALSE)
dir.create("./Analysis/Multiomics_Survival_Model/002_Multivariate_survival_analysis", 
           showWarnings = FALSE)

# Set parameters
Surv.cutoff.years = 20
Survival = "OS" # "DFS" or "OS"
version = "Micro_Low_ICR_High_vs_rest" # "Microbiome_HL_risk_in_ICR_High" # "Microbiome_HL_risk_in_ICR_Low"
                                            # Microbiome_HL_risk_in_ICR_Medium"

#"ICR_HML" # "Micr_ICR_High_Micro_ICRML"
# "Microbiome_ICR_combined_HML" # "Micro_Low_ICR_High_vs_rest" # "Micr_ICR_High_vs_rest" # "Microbiome_HL_risk_in_ICR_Medium_Low"
subset = "246" # "246" or "42"

if(Survival == "OS"){
  Status = "OS.Status"
  Time = "OS.Time"
}

if(Survival == "DFS"){
  Status = "DFS.Status"
  Time = "DFS.Time"
}

# Load data
load("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
Merged_dataset = merge(df, clinical_data, by = "Patient_ID")

Merged_dataset$var02_ICR_cluster = factor(Merged_dataset$var02_ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))

Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Alive")] = "0"
Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Dead")] = "1"

Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Disease Free")] = "0"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Event")] = "1"

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset[, Status] == "0", ]
TS.Alive$Status = TS.Alive[, Status]
TS.Alive$Time = TS.Alive[, Time]
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset[, Status] == "1", ]
TS.Dead$Status = TS.Dead[, Status]
TS.Dead$Time = TS.Dead[, Time]
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead[, Status][TS.Dead$Time> Y] = "0"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

if(subset == "246"){
  TS.Surv = TS.Surv[which(TS.Surv$Microbiome_Paired_TN == "yes"),]
}
if(subset == "288"){
  TS.Surv = TS.Surv[which(TS.Surv$Microbiome_T== "yes"),]
}
if(subset == "42"){
  TS.Surv = TS.Surv[which(TS.Surv$Microbiome_T == "yes"),]
  TS.Surv = TS.Surv[-which(TS.Surv$Microbiome_Paired_TN == "yes"),]
}

if(version == "Micro_Low_ICR_High_vs_rest"){
  TS.Surv$Combined_group = paste(TS.Surv$var02_ICR_cluster, TS.Surv$var16_Microbiome_risk_group, sep = "- ")
  TS.Surv$Combined_group[-which(TS.Surv$Combined_group == "ICR High- Low Risk")] = "Rest"
  TS.Surv$Combined_group = factor(TS.Surv$Combined_group, levels = c("ICR High- Low Risk", "Rest"))
  palette = c("red", "black")
  linetypes = c("solid", "solid")
  
  height = 5
  width = 6
}

if(version == "Microbiome_HL_risk_in_ICR_High"){
  TS.Surv$Combined_group = paste(TS.Surv$var02_ICR_cluster, TS.Surv$var16_Microbiome_risk_group, sep = "- ")
  TS.Surv = TS.Surv[which(TS.Surv$Combined_group %in% c("ICR High- Low Risk", "ICR High- High Risk")),]
  TS.Surv$Combined_group = factor(TS.Surv$Combined_group, levels = c("ICR High- Low Risk", "ICR High- High Risk"))
  palette = c("red", "#F20000")
  linetypes = c("solid", "dashed")
  
  height = 5
  width = 6
}


Group.of.interest = "Combined_group"

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[, Group.of.interest], conf.type = "log-log")

summary(mfit)

ggsurv_plot = ggsurvplot(mfit,
                         data = TS.Surv,
                         xlim = c(0, 90),
                         break.time.by = 30,
                         censor = TRUE,
                         risk.table = TRUE,
                         tables.y.text.col = TRUE,
                         tables.y.text = FALSE,
                         tables.height = 0.3,
                         tables.theme = theme_cleantable(),
                         #tables.col = "strata",
                         risk.table.pos = "out",
                         legend = "none",
                         ylab = paste0(Survival, " probability"),
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.y = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         linetype = linetypes,
                         #pval = TRUE,
                         palette = palette
)

TS.Surv$Combined_group = factor(TS.Surv$Combined_group, levels = c("ICR High- High Risk", "ICR High- Low Risk"))
uni_variate = coxph(formula = Surv(Time, Status) ~ Combined_group, data = TS.Surv)
summary = summary(uni_variate)
summary
logrank_p = summary$logtest[3]

#ggsurv_plot$plot <- ggsurv_plot$plot+ 
 # ggplot2::annotate("text", 
  #                  x = 1, y = 0.2, # x and y coordinates of the text
   #                 label = paste0("logrank p value = ", signif(logrank_p, 3)), 
    #                size = 4.5, hjust = 0)

dir.create("./Figures/Multiomics_Survival_Model/003_Combined_KM_Multiple_Groups", showWarnings = FALSE)
#png(paste0("./Figures/Multiomics_Survival_Model/003_Combined_KM_Multiple_Groups/003_KM_", Survival, "_", version, "_in_",
 #          subset, ".png"),
  #  res=600, height = height, width=width,unit="in")
#print(ggsurv_plot)
#dev.off()

pdf(paste0("./Figures/Multiomics_Survival_Model/003_Combined_KM_Multiple_Groups/003_KM_", Survival, "_", version, "_in_",
           subset, ".pdf"),
    height=4.5, width=5.2, onefile = FALSE, family = "ArialMT")
print(ggsurv_plot)
dev.off()




