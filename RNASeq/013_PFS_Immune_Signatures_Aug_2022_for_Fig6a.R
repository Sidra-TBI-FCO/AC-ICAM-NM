


# 
# Resulting HR table can be used for forest plot in script 014

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
required.packages <- c("stringr", "survival")
ipak(required.packages)

# Set parameters
exclude_adjuvant = "" #"adjuvant_treated_excluded" or ""
Surv.cutoff.years = 20

# Load data
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load("./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/Immune_signatures_clusters_Roelands.Rdata")

## Survival Analysis with all continuous variable in seperate script:
# Exclude adjuvant treated
if(exclude_adjuvant == "adjuvant_treated_excluded"){
  clinical_data = clinical_data[which(clinical_data$Adjuvant_treatment == "No"),]
}

rownames(clinical_data) = clinical_data$Patient_ID
rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 3)

i=1
for (i in 1:ncol(immune_sig_df)){
  col = colnames(immune_sig_df)[i]
  immune_sig_df[, col] = (immune_sig_df[, col] - min(immune_sig_df[, col]))/(max(immune_sig_df[,col])-min(immune_sig_df[,col]))
}

clinical_data = merge(clinical_data, immune_sig_df, by = "row.names")

HR_table = data.frame(Signature = colnames(immune_sig_df), p_value_logrank = NA, p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)
i=1
for (i in 1:ncol(immune_sig_df)){
  Group.of.interest = colnames(immune_sig_df)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  TS.Alive = clinical_data[clinical_data$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$DFS.Status == "Event", c("DFS.Status", "DFS.Time", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Disease Free"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Event"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  p_value_logrank = summary$sctest[3]
  HR_table$p_value_logrank[which(HR_table$Signature == Group.of.interest)] = p_value_logrank
  HR_table$p_value[which(HR_table$Signature == Group.of.interest)] = p_value
  HR_table$CI_lower[which(HR_table$Signature == Group.of.interest)] = CI_lower
  HR_table$CI_upper[which(HR_table$Signature == Group.of.interest)] = CI_upper
  HR_table$HR[which(HR_table$Signature == Group.of.interest)] = HR
}

dir.create("./Analysis/Trimmed_p/013_Immune_signatures_survival", showWarnings = FALSE)

rownames(HR_table) = gsub(".*\\- ", "", HR_table$Signature)
HR_table = HR_table[df$Full_name,]

save(HR_table, file = paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_Aug_2022_PFS_HR_table_immune_signatures_",
                             exclude_adjuvant, "_cutoff_", Surv.cutoff.years,".Rdata"))

write.csv(HR_table, file = "./Analysis/Trimmed_p/013_Immune_signatures_survival/013_Aug_2022_PFS_HR_table_immune_signatures__cutoff_20.csv",
          row.names = TRUE)



