
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "survminer", "survival", "coin"))

# Set parameters
Surv.cutoff.years = 20  
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full" 
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
Tissue = "T"


# Load data
load(file= paste0("./Analysis/Microbiome/020_Categories_by_median/020_Categorized_", Type, "_abundance_by_median_", Tissue, "_", Rank,".Rdata"))
colnames(abundance_cat)
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

clinical_data = clinical_data[which(clinical_data$Patient_ID %in% rownames(abundance_cat)),]

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset$ICR_HML = factor(Merged_dataset$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]

results = data.frame(Name = colnames(abundance_cat), p_val = NA, HR = NA, CI_lower = NA, CI_upper = NA)

i = 371 # 277 (=Fusicatenibacter) # 371 (= Ruminococcus 1) # 372 (=Ruminococcus 2) # 419 Fusobacterium
for(i in 1:ncol(abundance_cat)){
  micr = colnames(abundance_cat)[i]
  # add microbiome variable
  Merged_dataset$Microbiome_var = abundance_cat[,micr][match(Merged_dataset$Patient_ID, rownames(abundance_cat))]
  
  
  # time / event object creation
  Y = Surv.cutoff.years * 365
  TS.Alive = Merged_dataset[Merged_dataset$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", "ICR_HML", 
                                                                      "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
  TS.Alive$Time[TS.Alive$Time > Y] = Y
  
  TS.Dead = Merged_dataset[Merged_dataset$DFS.Status == "Event", c("DFS.Status", "DFS.Time", "ICR_HML",
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", "MSI", "Microbiome_var")]
  colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "Microbiome_var")
  TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
  TS.Dead$Status[which(TS.Dead$Time> Y)] = "Disease Free"
  TS.Dead$Time[TS.Dead$Time > Y] = Y
  
  TS.Surv = rbind (TS.Dead,TS.Alive)
  TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
  TS.Surv$Status <- TS.Surv$Status == "Event"
  
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
  
  results$p_val[which(results$Name == micr)] = logrank_pval
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
  #ggsurv_plot$plot <- ggsurv_plot$plot + labs(
   # title    = gsub(".*\\D5__", "", micr),                     
  #  subtitle = "Relative abundance by median",
  #  caption = gsub(".*\\D_1__", "", micr)
  #)
  #ggsurv_plot <- ggpar(
   # ggsurv_plot,
  #  font.title    = c(12, "plain", "black"),         
  #  font.subtitle = c(12, "plain", "black"), 
   # font.caption  = c(5, "plain", "black")
  #)
  
  if(table(TS.Surv$Microbiome_var)[1] < 5 | table(TS.Surv$Microbiome_var)[2] < 5){next}
  if(ncol(abundance_cat) < 30 | logrank_pval < 0.1){
    dir.create("./Figures/Microbiome/021_KM_relative_abundance_groups", showWarnings = FALSE)
     png(paste0("./Figures/Microbiome/021_KM_relative_abundance_groups/KM_PFS", Type, "_abundance_in_", Tissue, "_tissue_", Rank,
               "_", micr, ".png"),
      res=600,height=3.8,width=3.8,unit="in")  # set filename
    print(ggsurv_plot)
    dev.off()
  }
  
}

dir.create("./Analysis/Microbiome/021_KM_by_median", showWarnings = FALSE)
save(results, file = paste0("./Analysis/Microbiome/021_KM_by_median/PFS_KM_relative_abundance_by_median_", Tissue ,".Rdata"))


