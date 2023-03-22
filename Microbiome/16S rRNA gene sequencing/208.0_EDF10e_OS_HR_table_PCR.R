
# OS samples PCR 

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx", "stringr"))

# set parameters
Surv.cutoff.years = 20
cohort = "overlapping" # all # overlapping
survival = "OS" # OS  # PFS

# load data
load("./Processed_Data/Survival Data/JSREP_NT_clinical_data.Rdata")
pcr = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/PCR_validation/Microbiome validation data all WGS cohort complete.xlsx")
samples.list = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/PCR_validation/list_of_samples_each_batch.xlsx")

# subset PCR 
pcr = pcr[!is.na(pcr$PCR),]

pcr$Cohort = str_trim(pcr$Cohort, "both")

if(cohort == "all"){
  pcr = pcr[which(pcr$Cohort == "training"),]
}

if(cohort == "overlapping"){
  pcr = pcr[!is.na(pcr$`WGS.Species.=.Bromii`),]
}


# subset clinical data
clinical_data = clinical_data[which(clinical_data$ID %in% pcr$Sample),]

# add annotation 
clinical_data$bromi = pcr$Annotation[match(clinical_data$ID, pcr$Sample)]
clinical_data$bromi = factor(clinical_data$bromi, levels = c("Negative", "Positive"))

clinical_data <- clinical_data %>% relocate(bromi, .before = Type)


table(clinical_data$bromi)

Merged_dataset = clinical_data

### OS
# time / event object creation
Y = Surv.cutoff.years * 365

if (survival == "OS") {
  TS.Alive = Merged_dataset[Merged_dataset$OS.Status == "Alive", c("OS.Status", "OS.Time", "bromi")]
  colnames(TS.Alive) = c("Status","Time", "Group")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Merged_dataset[Merged_dataset$OS.Status == "Dead", c("OS.Status", "OS.Time", "bromi")]
  colnames(TS.Dead) = c("Status","Time", "Group")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Dead"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1) # remove patients with less then 1 day follow up time
}

if (survival == "PFS") {
  TS.Alive = Merged_dataset[Merged_dataset$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", "bromi")]
  colnames(TS.Alive) = c("Status","Time", "Group")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Merged_dataset[Merged_dataset$DFS.Status == "Event", c("DFS.Status", "DFS.Time", "bromi")]
  colnames(TS.Dead) = c("Status","Time", "Group")                     
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Disease Free"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Event"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)  # remove patients with less then 1 day follow up time
  
}

table(TS.Surv$Group)

fit = survfit(Surv(TS.Surv$Time/30.4, TS.Surv$Status) ~ TS.Surv[,"Group"], data = TS.Surv, conf.type = "log-log")

# HR 
uni_variate = coxph(formula = Surv(Time, Status) ~ Group, data = TS.Surv)
summary = summary(uni_variate)
HR = summary$conf.int[1]
CI_lower = summary$conf.int[3]
CI_upper = summary$conf.int[4]
p_value = summary$coefficients[5]
p_value_logrank = summary$logtest[3]
summary

results_df = data.frame(Survival = "PCR-Ruminococcus bromii", HR = NA, p_val = NA, CI_lower = NA, CI_upper = NA, p_val_logrank = NA)

results_df$HR = HR
results_df$p_val = p_value
results_df$CI_lower = CI_lower
results_df$CI_upper = CI_upper
results_df$p_val_logrank = p_value_logrank 

save(results_df, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/021_HR_tables_Median_and_continuous/HR_PCR_",survival,"_",cohort,"_samples_categorical.Rdata"))
