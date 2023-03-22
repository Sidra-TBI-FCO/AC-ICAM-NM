

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "survminer", "survival", "coin"))

# Set parameters
Surv.cutoff.years = 20  
Survival = "DFS" # "DFS" or "OS"
MSI = "MSI-H" # or "MSI-H"

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters 
rownames(Species_WGS)
Species = "s__Fusobacterium_nucleatum"  # "s__Bacteroides_fragilis"  #"s__Fusobacterium_nucleatum"

# Getting MSS and MSI-H
include_patients = MANTIS$Patient_ID[which(MANTIS$MSI == MSI)]

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

# Analysis
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% substring(colnames(Species_WGS), 1, 3)),]

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset$ICR_HML = factor(Merged_dataset$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]

Merged_dataset$vital_status[which(Merged_dataset$vital_status == "Alive")] = "0"
Merged_dataset$vital_status[which(Merged_dataset$vital_status == "Dead")] = "1"

Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Disease Free")] = "0"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Event")] = "1"

if(MSI == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$Patient_ID %in% include_patients),]
}

results = data.frame(Name = Species, p_val = NA, HR = NA, CI_lower = NA, CI_upper = NA)

  micr = Species
  # add microbiome variable
  Merged_dataset$Microbiome_abundance = Species_WGS[Species,][match(Merged_dataset$Patient_ID, 
                                                              substring(colnames(Species_WGS), 1, 3))]
  
  Merged_dataset$Microbiome_var = "High"
  Merged_dataset$Microbiome_var[which(Merged_dataset$Microbiome_abundance == 0)] = "Low"
  
  
  table(Merged_dataset$Microbiome_var)
  

  # time / event object creation
  Y = Surv.cutoff.years * 365
  TS.Alive = Merged_dataset[Merged_dataset[, Status] == "0", c(Status, Time_no_event, "ICR_HML", 
                                                                      "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Merged_dataset[Merged_dataset[, Status] == "1", c(Status, Time_event, "ICR_HML",
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  
  if(length(unique(TS.Surv$Microbiome_var)) == 1){next}
  TS.Surv$Microbiome_var = factor(TS.Surv$Microbiome_var, levels = c("Low", "High"))
  
  # survival curve
  msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
  mfit = survfit(msurv~TS.Surv$Microbiome_var,conf.type = "log-log")
  
  # Calculations (Needs manual adaptation!)
  mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
  pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
  pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
  
  Cal.Surv = TS.Surv
  mHR = coxph(formula = msurv ~ Cal.Surv$Microbiome_var,data = Cal.Surv, subset = Cal.Surv$Microbiome_var %in% c("Low", "High"))
  
  logrank_test(Surv(Time, Status) ~ Microbiome_var, data = TS.Surv)
  test = summary(mHR)
  logrank_pval = test$logtest[3]
  
  p_val = test$coefficients[5]
  HR = test$coefficients[2]
  CI_lower = test$conf.int[3]
  CI_upper = test$conf.int[4]
  
  results$p_val[which(results$Name == micr)] = p_val
  results$HR[which(results$Name == micr)] = HR 
  results$CI_lower[which(results$Name == micr)] = CI_lower
  results$CI_upper[which(results$Name == micr)] = CI_upper
  
  # plots
  ggsurv_plot = ggsurvplot(mfit,
                           xlim = c(0, 90),
                           break.time.by = 30,
                           data = TS.Surv,
                           censor = TRUE,
                           risk.table = TRUE,
                           tables.y.text.col = FALSE,
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
                           pval = TRUE,
                           palette = c("darkred", "darkblue")
  ) 

  
    dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_07_Kaplan_Meier_OS_by_WGS_Species", showWarnings = FALSE)
    pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_07_Kaplan_Meier_OS_by_WGS_Species/EDF8m_KM_", MSI, "_",
               "Kaplan_Meier_", Species, "_", Survival ,"_by_zero_cutoff.pdf"),
        height=3.8,width=3.8)  # set filename
    print(ggsurv_plot)
    dev.off()

    
    
