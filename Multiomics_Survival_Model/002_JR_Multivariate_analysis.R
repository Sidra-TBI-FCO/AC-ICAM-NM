
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
Survival = "OS" # "OS" or "DFS"
subset = "246"

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

Merged_dataset$var02_ICR_cluster[which(Merged_dataset$var02_ICR_cluster == "ICR Medium")] = NA
Merged_dataset$var02_ICR_cluster = factor(Merged_dataset$var02_ICR_cluster, levels = c("ICR Low", "ICR High"))

Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Alive")] = "0"
Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Dead")] = "1"

Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Disease Free")] = "0"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Event")] = "1"

if(subset == "246"){
  Merged_dataset = Merged_dataset[which(Merged_dataset$Microbiome_Paired_TN == "yes"),]
}

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

table(TS.Surv$var09_GIE_category)

TS.Surv$var09_GIE_category = factor(TS.Surv$var09_GIE_category, levels = c("non GIE", "GIE"))
Group.of.interest = "var09_GIE_category"

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[, Group.of.interest], conf.type = "log-log")

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
                         #linetype = linetypes,
                         #pval = TRUE,
                         palette = c("darkblue", "orange")
)

png(paste0("./Figures/Multiomics_Survival_Model/002_Multivariate_survival_analysis/", Group.of.interest,
           "_KM_", Survival, "_plot.png"),
    #res = 600, height = 6, width = 5.4, units = "in")
    #res = 600, height = 4.2, width = 5.4, units = "in")
    res=600, height = 3.8, width=4.2,unit="in")
print(ggsurv_plot)
dev.off()


# Univariate all
variables = grep("var", colnames(df), value = TRUE)

results = data.frame(variable = variables, logrank_p = NA, 
                     coxph_HR = NA, coxph_p = NA, coxph_CI_lower = NA,coxph_CI_upper = NA)
i=4
for (i in 1:length(variables)){
  variable_of_interest = variables[i]
  
  if(sum(is.na(TS.Surv[, variable_of_interest])) > 0){
    subset_surv = TS.Surv[-which(is.na(TS.Surv[, variable_of_interest])),]
  }else{
    subset_surv = TS.Surv
  }
  
  uni_variate = coxph(formula = Surv(Time, Status) ~ get(variable_of_interest), data = subset_surv)
  summary = summary(uni_variate)
  logrank_p = summary$logtest[3]
  HR = summary$coefficients[2]
  coxph_p = summary$coefficients[5]
  CI_lower = summary$conf.int[3]
  CI_upper = summary$conf.int[4]
  
  results$logrank_p[i] = signif(logrank_p, 3)
  results$coxph_HR[i] = signif(HR, 3)
  results$coxph_CI_lower[i] = signif(CI_lower, 3)
  results$coxph_CI_upper[i] = signif(CI_upper, 3)
  results$coxph_p[i] = signif(coxph_p, 3)
  
}

results = results[order(results$coxph_p, decreasing = FALSE),]

write.csv(results, file = paste0("./Analysis/Multiomics_Survival_Model/002_Multivariate_survival_analysis/Univariate_coxph_all_variables_", 
                                 subset, "_", Survival, ".csv"),
         row.names = FALSE)

TS.Surv$varnew_ICR_H_vs_ML = factor(TS.Surv$var20_ICR_ordinal)
levels(TS.Surv$varnew_ICR_H_vs_ML) = c("ICR M-L", "ICR M-L", "ICR High")

TS.Surv$varnew2_Stage4_vs_Rest = factor(TS.Surv$var17_Stage_ordinal)
levels(TS.Surv$varnew2_Stage4_vs_Rest) = c("Stage I-III", "Stage I-III", "Stage I-III", "Stage IV")

TS.Surv = TS.Surv[which(TS.Surv$Microbiome_Paired_TN == "yes"),]

colnames(TS.Surv)

# "var01_ICRscore"                                           "var02_ICR_cluster"                                    "var03_ICR_UP_DOWN"                                   
#[10] "var04_MiXCR_clonality"                                "var05_MiXCR_clonality_group"                          "var06_Adaptive_TCR_clonality"                        
#[13] "var07_Adaptive_TCR_clonality_group"                   "var08_GIE_value"                                      "var09_GIE_category"                                  
#[16] "var10_IES"                                            "var11_Nonsynonymous_Mutation_count"                   "var12_Hypermutation_status"                          
#[19] "var13_Neoantigen_count"                               "var14_Ratio_Neoantigen_count_Nonsynonymous_Mut_count" "var15_Microbiome_risk_score"                         
#[22] "var16_Microbiome_risk_group"                          "var17_Stage_ordinal"                                  "var18_Stage_categorical"                             
#[25] "var19_MSI_MANTIS"                                     "var20_ICR_ordinal"
# "var21_ICR_H_vs_ML"                                    "var22_gender", "var23_CMS",  "var22_Stage4_vs_Rest" "age_at_initial_pathologic_diagnosis"


# Multivariate
multi_variate = coxph(formula = Surv(Time, Status) ~ var20_ICR_ordinal + var16_Microbiome_risk_group + var09_GIE_category + var17_Stage_ordinal  + 
                        var21_age + var19_MSI_MANTIS + var23_CMS, data = TS.Surv)
summary = summary(multi_variate)
summary



multi_variate = coxph(formula = Surv(Time, Status) ~ var20_ICR_ordinal + var16_Microbiome_risk_group + var09_GIE_category + var17_Stage_ordinal  + 
                        var21_age, data = TS.Surv)
summary = summary(multi_variate)
summary




multi_variate = coxph(formula = Surv(Time, Status) ~ var03_ICR_UP_DOWN + var16_Microbiome_risk_group + var09_GIE_category + var18_Stage_categorical, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var21_ICR_H_vs_ML + var16_Microbiome_risk_group + var09_GIE_category + var18_Stage_categorical, data = TS.Surv)


multi_variate = coxph(formula = Surv(Time, Status) ~ var20_ICR_ordinal + var16_Microbiome_risk_group + var18_Stage_categorical + var09_GIE_category, data = TS.Surv)
multi_variate = coxph(formula = Surv(Time, Status) ~ var20_ICR_ordinal + var16_Microbiome_risk_group + var09_GIE_category, data = TS.Surv)
multi_variate = coxph(formula = Surv(Time, Status) ~ var20_ICR_ordinal + var16_Microbiome_risk_group + var17_Stage_ordinal + 
                        var09_GIE_category + age_at_initial_pathologic_diagnosis, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var21_ICR_H_vs_ML + var16_Microbiome_risk_group + var09_GIE_category + var17_Stage_ordinal, data = TS.Surv)

#TS.Surv = TS.Surv[which(TS.Surv$Microbiome_Paired_TN == "yes"),] 
multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + var21_ICR_H_vs_ML + var09_GIE_category + var18_Stage_categorical, data = TS.Surv)
multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + var21_ICR_H_vs_ML + var09_GIE_category + var17_Stage_ordinal, data = TS.Surv)
summary = summary(multi_variate)
summary 

multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + var21_ICR_H_vs_ML, data = TS.Surv)
summary = summary(multi_variate)
summary 

multi_variate = coxph(formula = Surv(Time, Status) ~ var21_ICR_H_vs_ML + var16_Microbiome_risk_group + var09_GIE_category +
                        var18_Stage_categorical, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var09_GIE_category + var16_Microbiome_risk_group, data = TS.Surv)
summary = summary(multi_variate)
summary 

multi_variate = coxph(formula = Surv(Time, Status) ~ var09_GIE_category + var16_Microbiome_risk_group + var21_ICR_H_vs_ML + var22_Stage4_vs_Rest, data = TS.Surv)
summary = summary(multi_variate)
summary 

multi_variate = coxph(formula = Surv(Time, Status) ~  var09_GIE_category + 
                        var03_ICR_UP_DOWN, data = TS.Surv)


multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + var09_GIE_category +
                        var03_ICR_UP_DOWN, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + var09_GIE_category +
                        var02_ICR_cluster, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var16_Microbiome_risk_group + 
                        var02_ICR_cluster, data = TS.Surv)

multi_variate = coxph(formula = Surv(Time, Status) ~ var15_Microbiome_risk_score +  var09_GIE_category +
                        var01_ICRscore, data = TS.Surv)

summary = summary(multi_variate)
summary 

TS.Surv$var03_ICR_UP_DOWN
TS.Surv$var09_GIE_category
TS.Surv$var16_Microbiome_risk_group
TS.Surv$var17_Stage_ordinal

