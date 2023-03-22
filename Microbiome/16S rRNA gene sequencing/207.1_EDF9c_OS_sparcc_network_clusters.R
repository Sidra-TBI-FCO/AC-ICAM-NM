
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival"))

# Set parameters
Surv.cutoff.years = 20
Survival = "OS" # "DFS" or "OS"
matrix = "RA" # RA # OTU
version = "without" # without # subtract

# load data
load("./Processed_Data/Survival Data/JSREP_NT_clinical_data.Rdata")
load(paste0("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/network/network_clusters_df_",version,"_Prevotella2_",matrix,".Rdata"))

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

clinical_data = clinical_data[which(clinical_data$ID %in% df$paitent_ID),]

colnames(df)[1] = "ID"

clinical_data = merge(clinical_data, df, "ID")

rownames(df) = df$ID
df$ID = NULL

clinical_data$ID == rownames(df)

clinical_data$vital_status[which(clinical_data$vital_status == "Alive")] = "0"
clinical_data$vital_status[which(clinical_data$vital_status == "Dead")] = "1"

clinical_data$DFS.Status[which(clinical_data$DFS.Status == "Disease Free")] = "0"
clinical_data$DFS.Status[which(clinical_data$DFS.Status == "Event")] = "1"

results = data.frame(Cluster = colnames(df), HR = NA, CI_lower = NA, CI_upper = NA, 
                     cox_pval = NA, logrank_pval = NA)

i = 4
for (i in 1:ncol(df)) {
  clus = df[,i]
  cluster.df = data.frame(ID = rownames(df), cluster = clus, vital_status = clinical_data$vital_status, last_contact_days_to = clinical_data$last_contact_days_to, 
                          death_days_to = clinical_data$death_days_to, DFS.Status = clinical_data$DFS.Status, DFS.Time = clinical_data$DFS.Time)
  
  m.cluster = median(cluster.df$cluster)
  
  cluster.df$group = NA
  cluster.df$group[which(cluster.df$cluster <= m.cluster)] = "Low abundance"
  cluster.df$group[which(cluster.df$cluster > m.cluster)] = "High abundance"
  cluster.df$group = factor(cluster.df$group, levels = c("Low abundance", "High abundance"))
  
  table(cluster.df$group)
  
  # time / event object creation
  Y = Surv.cutoff.years * 365
  
  TS.Alive =  cluster.df[cluster.df[, Status] == "0", c(Status, Time_no_event, "group")]
  colnames(TS.Alive) = c("Status","Time", "Group")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = cluster.df[cluster.df[, Status] == "1", c(Status, Time_event, "group")]
  colnames(TS.Dead) = c("Status","Time", "Group")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "1"
  TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time
  
  mfit = survfit(Surv(TS.Surv$Time/30.4, TS.Surv$Status) ~ TS.Surv[,"Group"], data = TS.Surv, conf.type = "log-log")
  
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
                           palette = c("darkblue", "darkred")
  ) 
  
  
  #dir.create("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/OS", showWarnings = FALSE)
  svg(paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/OS/",Survival,"_cluster_",i,"_",version,"_Prevotella2_",matrix,".svg"),
      height = 3, width = 5) # set filename
  print(ggsurv_plot)
  dev.off()
  
  # HR 
  uni_variate = coxph(formula = Surv(Time, Status) ~ Group, data = TS.Surv)
  summary = summary(uni_variate)
  HR = summary$conf.int[1]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  p_value = summary$coefficients[5]
  p_value_logrank = summary$logtest[3]
  summary
  
  
  results$HR[which(results$Cluster == colnames(df)[i])] = round(HR, 3)
  results$CI_lower[which(results$Cluster == colnames(df)[i])] = round(CI_lower, 3)
  results$CI_upper[which(results$Cluster == colnames(df)[i])] = round(CI_upper, 3)
  results$cox_pval[which(results$Cluster == colnames(df)[i])] = p_value
  results$logrank_pval[which(results$Cluster == colnames(df)[i])] = p_value_logrank
  
}

#dir.create("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/OS", showWarnings = FALSE)
save(results, file = paste0("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/OS/sparcc_HR_table_",Survival, "_network_clusters_",version,"_Prevotella2_",matrix,".Rdata"))

