


# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "stringr", "survminer", "survival", "coin", "survival", "openxlsx"))

# Set parameters
Surv.cutoff.years = 20
Dataset = "246" # "288" or "246" or "42" or "PCR subset"
IES = "" # "IES4" or "IES1_2_3"
Survival = "OS" # "OS" or "DFS"
risk_label = "risk_lab" # "risk_lab" or "risk_without_ruminococcus_lab" or "risk_genus_subset_lab" 

# 
type_label = gsub("risk", "prediction", risk_label)
continuous_risk_score = gsub("_lab", "", type_label) #prediction_without_ruminococcus" # "prediction" or "prediction_without_ruminococcus" or "prediction_subset_genus"

if(Survival == "OS"){
  Status = "vital_status"
  Time_event = "death_days_to"
  Time_no_event = "last_contact_days_to"
}

if(Survival == "DFS"){
  Status = "DFS.Status"
  Time_event = "DFS.Time"
  Time_no_event = "DFS.Time"
}


# Load data
load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
train_risk_scores_micr = Training
valid_risk_scores_micr = Validation

if(Dataset == "288"){
  all_risk_scores_micr = rbind(train_risk_scores_micr, valid_risk_scores_micr)
}
if(Dataset == "246"){
  all_risk_scores_micr = rbind(train_risk_scores_micr)
}
if(Dataset == "42"){
  all_risk_scores_micr = valid_risk_scores_micr
}

if(Dataset == "PCR subset"){
  pcr = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/Copy of AC-ICAM_Ruminococcus1_and_2__246_samples_relative_abundance_values_Chris7-8-22.xlsx")
  pcr = pcr[!is.na(pcr$BROMII.PCR.results),]
  pcr_Patient_ID = substring(pcr$Sample, 1, 3)
  all_risk_scores_micr = train_risk_scores_micr[which(train_risk_scores_micr$Patient_ID %in% pcr_Patient_ID),]
}

colnames(all_risk_scores_micr)[which(colnames(all_risk_scores_micr) == "prediction_subset_genus")] = "prediction_genus_subset"

load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset_include_medium_cat.Rdata")
GIE_df = Merged_dataset

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% all_risk_scores_micr$Patient_ID),]

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset$ICR_HML = factor(Merged_dataset$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]
Merged_dataset$IES = IES_df$IES[match(Merged_dataset$Patient_ID, IES_df$Patient_ID)]

Merged_dataset$GIE = GIE_df$Immunoediting_score[match(Merged_dataset$Patient_ID, GIE_df$Patient_ID)]
Merged_dataset$GIE_cat = GIE_df$Immunoedited[match(Merged_dataset$Patient_ID, GIE_df$Patient_ID)]

if(IES == "IES4"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$IES == "IES4"),]
}
if(IES == "IES1_2_3"){
  Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$IES)),]
  Merged_dataset = Merged_dataset[-which(Merged_dataset$IES == "IES4"),]
}
#Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$IES)),]

Merged_dataset$vital_status[which(Merged_dataset$vital_status == "Alive")] = "0"
Merged_dataset$vital_status[which(Merged_dataset$vital_status == "Dead")] = "1"

Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Disease Free")] = "0"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Event")] = "1"

Merged_dataset = merge(Merged_dataset, all_risk_scores_micr, by = "Patient_ID")

Merged_dataset$Risk_score_micr = all_risk_scores_micr[,risk_label][match(Merged_dataset$Patient_ID,
                                                                         all_risk_scores_micr$Patient_ID)]

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset[, Status] == "0", c(Status, Time_no_event, "ICR_HML", 
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", 
                                                                    "Risk_score_micr", "IES", "GIE", "GIE_cat", continuous_risk_score)]
colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Risk_score_micr", 
                       "IES", "GIE", "GIE_cat", "continuous_risk_score")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset[, Status] == "1", c(Status, Time_event, "ICR_HML",
                                                                  "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", 
                                                                  "Risk_score_micr", "IES", "GIE", "GIE_cat", continuous_risk_score)]
colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Risk_score_micr", "IES", "GIE", "GIE_cat", "continuous_risk_score")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"

if(length(unique(TS.Surv$Microbiome_var)) == 1){next}

table(TS.Surv$IES)
TS.Surv$IES = factor(TS.Surv$IES, levels = c("IES1", "IES2", "IES3", "IES4"))

TS.Surv$GIE_cat = factor(TS.Surv$GIE_cat, levels = c("less immunoedited", "immunoedited"))
table(TS.Surv$Risk_score_micr)

TS.Surv$Combined_group = paste(TS.Surv$ICR_cluster, TS.Surv$GIE_cat, TS.Surv$Risk_score_micr, sep = "+")
test = data.frame(table(TS.Surv$Combined_group))

dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier", showWarnings = FALSE)

palette = c("#e038d5", "#182ff5")
height = 3.8
width = 4.2

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv$Risk_score_micr,conf.type = "log-log")

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
                         ylab = paste0(Survival, " probability"),
                         xlab = "Time in months",
                         fontsize = 4.5,
                         font.x = 18,
                         font.y = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #linetype = linetypes,
                         #pval = TRUE,
                         palette = palette
)

TS.Surv$Risk_score_micr = factor(TS.Surv$Risk_score_micr, levels = c("Low Risk","High Risk"))

uni_variate = coxph(formula = Surv(Time, Status) ~ Risk_score_micr, data = TS.Surv)
summary = summary(uni_variate)
logp = signif(summary$logtest[3], 3)
HR = signif(summary$coefficients[2], 3)
CI_low = signif(summary$conf.int[3], 3)
CI_high = signif(summary$conf.int[4], 3)
p_cox = signif(summary$coefficients[5], 3)
continuous = coxph(formula = Surv(Time, Status) ~ continuous_risk_score, data = TS.Surv)
cont_summary = summary(continuous)
cont_pvalue = signif(cont_summary$coefficients[5], 3)

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/101_KM_group_Micr/101_", 
           Survival, "_KM_group_", risk_label, "_", Dataset, "_", IES, ".pdf"),
    height = height, width=width, onefile = FALSE)
print(ggsurv_plot)
dev.off()

summary

ggsurv_plot$plot <- ggsurv_plot$plot+ 
  ggplot2::annotate("text", 
                    x = 1, y = 0.2, # x and y coordinates of the text
                    label = paste0(logp, "\nHR = ", HR, " (", CI_low, "-", CI_high, "), p = ", p_cox,
                                   "\nContinuous p= ", cont_pvalue), 
                    size = 4.5, hjust = 0)

dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/101_KM_group_Micr", showWarnings = FALSE)

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/101_KM_group_Micr/101_", 
          Survival, "_KM_group_", risk_label, "_", Dataset, "_", IES, ".png"),
    res=600, height = height, width=width,unit="in")
print(ggsurv_plot)
dev.off()


table(TS.Surv$Risk_score_micr)

TS.Surv$Risk_score_micr = factor(TS.Surv$Risk_score_micr, levels = c("High Risk", "Low Risk"))

uni_variate = coxph(formula = Surv(Time, Status) ~ Risk_score_micr, data = TS.Surv)
summary = summary(uni_variate)
logp = summary$logtest[3]
HR = summary$coefficients[2]
CI_low = summary$conf.int[3]
CI_high = summary$conf.int[4]
p_cox = summary$coefficients[5]
continuous = coxph(formula = Surv(Time, Status) ~ continuous_risk_score, data = TS.Surv)
cont_summary = summary(continuous)
cont_pvalue =cont_summary$coefficients[5]
summary

if(Dataset == "246"){
  Name = "AC-ICAM - Tumor"
}

if(Dataset == "42"){
  Name = "ICAM42 - Tumor"
}

results_df = data.frame(Survival = Name, HR = HR, p_val = p_cox, 
                        CI_lower = CI_low, CI_upper = CI_high, p_val_logrank = logp)

dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier",
           showWarnings = FALSE)

save(results_df, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/", 
     "101_Kaplan_Meier/HR_table_", Survival, "_", Dataset, ".Rdata"))
