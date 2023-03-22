


# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize",
                                   "dendsort", "stringr")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
Surv.cutoff.years = 20
geneset = "ConsensusTME_COAD"

# Load
clinical_data = read.csv("./Processed_Data/External_Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                         stringsAsFactors = FALSE) # Read in the clinical data file
load(paste0("./Analysis/Deconvolution_and_GSEA/", geneset,"_ES.Rdata"))
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")

# Prepare data

rownames(clinical_data) = clinical_data$bcr_patient_barcode

immune_sig_df = data.frame(t(ES))

# Add ICR score
table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(ES)),]
immune_sig_df$ICR_score = table_cluster_assignment$ICRscore[match(rownames(immune_sig_df),
                                                                  rownames(table_cluster_assignment))]

rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 12)

i=1
for (i in 1:ncol(immune_sig_df)){
  col = colnames(immune_sig_df)[i]
  immune_sig_df[, col] = (immune_sig_df[, col] - min(immune_sig_df[, col]))/(max(immune_sig_df[,col])-min(immune_sig_df[,col]))
}

clinical_data = merge(clinical_data, immune_sig_df, by = "row.names")

HR_table = data.frame(Signature = colnames(immune_sig_df), p_value_logrank = NA, p_value = NA, HR = NA, CI_lower = NA, CI_upper = NA)
i=2
for (i in 1:ncol(immune_sig_df)){
  Group.of.interest = colnames(immune_sig_df)[i]
  Y = Surv.cutoff.years * 365
  # time / event object creation
  TS.Alive = clinical_data[clinical_data$PFI == "0", c("PFI", "PFI.time", Group.of.interest)]
  colnames(TS.Alive) = c("Status","Time", Group.of.interest)
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = clinical_data[clinical_data$PFI == "1", c("PFI", "PFI.time", Group.of.interest)]
  colnames(TS.Dead) = c("Status","Time", Group.of.interest)
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
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

dir.create("./Analysis/013_Immune_signatures_survival", showWarnings = FALSE)
save(HR_table, file = paste0("./Analysis/013_Immune_signatures_survival/013_", geneset, "_PFI_HR_table_",
                             "cutoff_", Surv.cutoff.years,".Rdata"))

