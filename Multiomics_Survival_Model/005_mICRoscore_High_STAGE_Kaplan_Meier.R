

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

ipak(c("forestplot", "survminer", "survival"))

# Set parameters
Group.of.interest = "ajcc_pathologic_tumor_stage"
x = 2
Stages = "all stages" #"all stages" #"all stages" #"all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name =  "All" #"Stage III" # "All" #"Stage III&IV" #"All"
Surv.cutoff.years = 20
Survival = "DFS" # "DFS" or "OS"
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

# Overview subgroups
table(Merged_dataset$var02_ICR_cluster, Merged_dataset$var16_Microbiome_risk_group, exclude = NULL)
Merged_dataset$Combined_group = paste(Merged_dataset$var02_ICR_cluster, Merged_dataset$var16_Microbiome_risk_group, sep = "- ")

Merged_dataset = Merged_dataset[which(Merged_dataset$Combined_group == "ICR High- Low Risk"),]

table(Merged_dataset$ajcc_pathologic_tumor_stage)
if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}
table(Merged_dataset$ajcc_pathologic_tumor_stage)

### KM
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
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1) 

if(subset == "246"){
  TS.Surv = TS.Surv[which(TS.Surv$Microbiome_Paired_TN == "yes"),]
}

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv$ajcc_pathologic_tumor_stage,conf.type = "log-log")

mfit_sum = summary(mfit)
mfit_sum

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

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
                         ylab = "",
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #linetype = linetypes,
                         #pval = TRUE,
                         palette = c("purple", "violet", "blue", "red")
)

dir.create("./Figures/Multiomics_Survival_Model/005_mICRoScore_High_Stages_Kaplan_Meiers", showWarnings = FALSE)
png(paste0("./Figures/Multiomics_Survival_Model/005_mICRoScore_High_Stages_Kaplan_Meiers/005_",
           Survival, "_", Group.of.interest, "_", subset, "_in_mICRoScore_High.png"),
    #res = 600, height = 5, width = 6, units = "in")
    #res = 600, height = 4.2, width = 5.4, units = "in")
    res=600, height = 3.8, width=4.2,unit="in")

print(ggsurv_plot)
dev.off()


