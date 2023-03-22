

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "stringr", "survminer", "survival", "coin", "survival"))

# Set parameters
Surv.cutoff.years = 20
Dataset = "TCGA+42_Validation" # "TCGA+42_Validation"

# Load data
load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
train_risk_scores_micr = Training
valid_risk_scores_micr = Validation
load("./Processed_Data/Microbiome/External/9_August/TCGA_Microbiome_Risk_Scores.Rdata")

TCGA_valid_risk_score_micr = TCGA_microbiome_risk

train_risk_scores_micr$patientid = str_pad(train_risk_scores_micr$samples, 3, pad = "0")
valid_risk_scores_micr$patientid = str_pad(valid_risk_scores_micr$samples, 3, pad = "0")

TCGA_valid_risk_score_micr = TCGA_valid_risk_score_micr[, c("samples", "prediction")]
colnames(TCGA_valid_risk_score_micr) = c("patientid", "risk_score")

valid_risk_scores_micr = valid_risk_scores_micr[, c("patientid", "prediction")]
colnames(valid_risk_scores_micr) = c("patientid", "risk_score")

if(Dataset == "TCGA+42_Validation"){
  all_risk_scores_micr = rbind(valid_risk_scores_micr, TCGA_valid_risk_score_micr)
}

load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

clinical_data_TCGA = read.csv("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/External_Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                              stringsAsFactors = FALSE)
colnames(clinical_data_TCGA)[which(colnames(clinical_data_TCGA) == "bcr_patient_barcode")] = "Patient_ID"

clinical_data_TCGA = clinical_data_TCGA[, c("Patient_ID", "OS", "OS.time", "PFI", "PFI.time")]
clinical_data = clinical_data[, c("Patient_ID", "OS.Status", "OS.Time", "DFS.Status", "DFS.Time")]
colnames(clinical_data) = c("Patient_ID", "OS.Status", "OS.Time", "PFI.Status", "PFI.Time")
colnames(clinical_data_TCGA) = c("Patient_ID", "OS.Status", "OS.Time", "PFI.Status", "PFI.Time")

# Merge the TCGA and AC-ICAM clinical data
Merged_dataset = rbind(clinical_data_TCGA, clinical_data)

Merged_dataset$Risk_score_micr = all_risk_scores_micr$risk_score[match(Merged_dataset$Patient_ID,
                                                                       all_risk_scores_micr$patientid)]

Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$Risk_score_micr)),]

# Checks and harmonization
table(Merged_dataset$OS.Status)
Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == 1)] = "Dead"
Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == 0)] = "Alive"
table(Merged_dataset$OS.Status) # Alive  Dead 
#   102    57

table(Merged_dataset$PFI.Status)
Merged_dataset$PFI.Status[which(Merged_dataset$PFI.Status == 1)] = "Event"
Merged_dataset$PFI.Status[which(Merged_dataset$PFI.Status == 0)] = "Disease Free"
table(Merged_dataset$PFI.Status)  # Disease Free   Event
# 103            56

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$OS.Status == "Alive", c("OS.Status", "OS.Time", "Risk_score_micr")]
colnames(TS.Alive) = c("Status","Time", "Risk_score_micr")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$OS.Status == "Dead", c("OS.Status", "OS.Time", "Risk_score_micr")]
colnames(TS.Dead) = c("Status","Time", "Risk_score_micr")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"

TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

TS.Surv$Risk_score_micr_cat = NA
TS.Surv$Risk_score_micr_cat[which(TS.Surv$Risk_score_micr < 0)] = "Microbiome Low Risk"
TS.Surv$Risk_score_micr_cat[which(TS.Surv$Risk_score_micr >= 0)] = "Microbiome High Risk"
table(TS.Surv$Risk_score_micr_cat)

dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/027_Kaplan_Meier", showWarnings = FALSE)

palette = c("#e038d5", "#182ff5")
height = 3.8
width = 4.2

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv$Risk_score_micr_cat,conf.type = "log-log")

ggsurv_plot = ggsurvplot(mfit,
                         data = TS.Surv,
                         xlim = c(0, 90),
                         break.time.by = 30,
                         censor = TRUE,
                         risk.table = TRUE,
                         tables.y.text.col = TRUE,
                         tables.y.text = TRUE,
                         tables.height = 0.3,
                         tables.theme = theme_cleantable(),
                         #tables.col = "strata",
                         risk.table.pos = "out",
                         legend = "none",
                         risk.table.y.text = FALSE,
                         ylab = "",
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #linetype = linetypes,
                         #pval = TRUE,
                         palette = palette
)

dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/027_KM_group_Micr", showWarnings = FALSE)

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/027_KM_group_Micr/New_027_KM_group_Micr_", Dataset,".png"),
    res=600, height = height, width=width,unit="in")
print(ggsurv_plot)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/027_KM_group_Micr/New_027_KM_group_Micr_", Dataset,".pdf"),
    height = height, width=width)
print(ggsurv_plot)
dev.off()

table(TS.Surv$Risk_score_micr_cat)

#TS.Surv$Risk_score_micr_cat = factor(TS.Surv$Risk_score_micr_cat, levels = c("Microbiome Low Risk",
 #                                                                            "Microbiome High Risk"))

#uni_variate = coxph(formula = Surv(Time, Status) ~ Risk_score_micr_cat, data = TS.Surv)
#summary(uni_variate)
#logrank_test(Surv(Time, Status) ~ Risk_score_micr_cat, data = TS.Surv)

TS.Surv$Risk_score_micr = TS.Surv$Risk_score_micr_cat
TS.Surv$Risk_score_micr = factor(TS.Surv$Risk_score_micr, levels = c("Microbiome High Risk", "Microbiome Low Risk"))

uni_variate = coxph(formula = Surv(Time, Status) ~ Risk_score_micr, data = TS.Surv)
summary = summary(uni_variate)
logp = summary$logtest[3]
HR = summary$coefficients[2]
CI_low = summary$conf.int[3]
CI_high = summary$conf.int[4]
p_cox = summary$coefficients[5]
#continuous = coxph(formula = Surv(Time, Status) ~ continuous_risk_score, data = TS.Surv)
#cont_summary = summary(continuous)
#cont_pvalue =cont_summary$coefficients[5]
summary

if(Dataset == "TCGA+42_Validation"){
  Name = "ICAM42_TCGA_COAD - Tumor"
}


results_df = data.frame(Survival = Name, HR = HR, p_val = p_cox, 
                        CI_lower = CI_low, CI_upper = CI_high, p_val_logrank = logp)

dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier",
           showWarnings = FALSE)

save(results_df, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/", 
                               "101_Kaplan_Meier/HR_table_OS_", Dataset, ".Rdata"))

