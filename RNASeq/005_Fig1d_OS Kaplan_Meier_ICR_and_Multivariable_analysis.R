
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr", "survminer", "coin")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 20                                 # SET cut-off
Group.of.interest = "ICR_HML"                          # "tumour_anatomic_site" or "ICR_HML"
Stages = "all stages" #"all stages" # "all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name = "All" #"All" #"Stage I&II" #"All" # "Stage III&IV"
CMS = "CMS4" # "" for all, or "CMS4"
exclude = c("Conpair_lower_90_percent", "non-epithelial")  #c("Conpair_lower_90_percent", "non-epithelial") #"Conpair_lower_90_percent" 
                                                          # c("Conpair_lower_90_percent", "non-epithelial")

# Create folders
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Logfiles/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p/Survival Analysis", showWarnings = FALSE)

# Load data
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
immunoediting_df = Merged_dataset
rm(list = "Merged_dataset")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Add ICR as a variable and assign ICR cluster according to table cluster assignment
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("ICR Low", "ICR Medium", "ICR High"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Merged_dataset$MSI = MANTIS$MSI[match(Merged_dataset$Patient_ID, MANTIS$Patient_ID)]

Merged_dataset$GIE_value = immunoediting_df$Immunoediting_score[match(Merged_dataset$Patient_ID,
                                                                      immunoediting_df$Patient_ID)]

Merged_dataset$GIE_group = immunoediting_df$Immunoedited[match(Merged_dataset$Patient_ID,
                                                               immunoediting_df$Patient_ID)]

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}

if(exclude == ""){
}else{
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}

if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}

table(Merged_dataset$ajcc_pathologic_tumor_stage)
table(Merged_dataset$ajcc_pathologic_tumor_stage)
table(Merged_dataset$ICR_HML)

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$vital_status == "Alive", c("vital_status", "last_contact_days_to", Group.of.interest, 
                                                                    "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", 
                                                                    "MSI", "GIE_value", "GIE_group")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "GIE_value", "GIE_group")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$vital_status == "Dead", c("vital_status", "death_days_to", Group.of.interest, 
                                                                  "ICRscore", "ajcc_pathologic_tumor_stage", "CMS", "age_at_initial_pathologic_diagnosis", 
                                                                  "MSI", "GIE_value", "GIE_group")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "ICRscore", "pathologic_stage", "CMS", "Age", "MSI", "GIE_value", "GIE_group")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[,Group.of.interest],conf.type = "log-log")

mfit_sum = summary(mfit)
mfit_sum

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
#Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
#Cal.Surv[,Group.of.interest] = relevel(Cal.Surv[,Group.of.interest], "ICR3")
Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("ICR Low", "ICR Medium", "ICR High"))
mHR.extract = extract(mHR, include.aic = TRUE,
                      include.rsquared = TRUE, include.maxrs=TRUE,
                      include.events = TRUE, include.nobs = TRUE,
                      include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)

logrank_test(Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv[which(TS.Surv[, Group.of.interest] %in% c("ICR High", "ICR Medium", "ICR Low")),])
summary = summary(mHR)
summary

Cal.Surv$ICR_ordinal = Cal.Surv$ICR_HML
levels(Cal.Surv$ICR_ordinal) = c("1", "2", "3")
Cal.Surv$ICR_ordinal = as.character(Cal.Surv$ICR_ordinal)
Cal.Surv$ICR_ordinal = as.numeric(Cal.Surv$ICR_ordinal)
ordinal_coxph = coxph(formula = msurv ~ Cal.Surv$ICR_ordinal, data = Cal.Surv)
summary = summary(ordinal_coxph)
summary$logtest
summary


# plot
TS.Surv[,Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels = c("ICR High", "ICR Medium", "ICR Low"))
mfit = survfit(msurv~TS.Surv[,Group.of.interest],conf.type = "log-log")

ggsurv_plot = ggsurvplot(mfit,
                         xlim = c(0, 90),
                         break.time.by = 30,
                         data = TS.Surv,
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
                         #fontsize = 4.5,
                         font.main = 12,
                         font.x = 12,
                         font.y = 12,
                         font.tickslab = 12,
                         font.caption = 12,
                         font.legend = 12,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #pval = TRUE,
                         palette = c("red", "green", "blue")
)
  
dir.create("./Figures/Trimmed_p/005_ggsurv_plots", showWarnings = FALSE)
pdf(paste0("./Figures/Trimmed_p/005_ggsurv_plots/Fig1d_", CMS,"_Overall_Survival_", Group.of.interest,"_Surv_cutoff_years_",Surv.cutoff.years,".pdf"),
    height=4,width=4.8, onefile = FALSE, family = "ArialMT")  # set filename
print(ggsurv_plot)
dev.off()

TS.Surv$MSI = factor(TS.Surv$MSI, levels = c("MSS", "MSI-H"))
TS.Surv$ICR_ordinal = TS.Surv$ICR_HML
levels(TS.Surv$ICR_ordinal) = c("1", "2", "3")
TS.Surv$ICR_ordinal = as.character(TS.Surv$ICR_ordinal)
TS.Surv$ICR_ordinal = as.numeric(TS.Surv$ICR_ordinal)
TS.Surv$ICR_HML[which(TS.Surv$ICR_HML == "ICR Medium")] = NA

# Multivariate analysis (Supplementary Table 2)

if(cox == "cox1"){
  # Cox regression (continuous)
  uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
  summary(uni_variate_ICRscore)
  
  # Cox regression MSI
  uni_variate_MSI = coxph(formula = Surv(Time, Status) ~ MSI, data = TS.Surv)
  summary(uni_variate_MSI)
  
  # Univariate stage
  TS.Surv$pathologic_stage = as.numeric(TS.Surv$pathologic_stage)
  uni_variate_stage = coxph(formula = Surv(Time, Status) ~ pathologic_stage, data = TS.Surv)
  summary(uni_variate_stage)
  
  # Age
  TS.Surv$Age = as.numeric(TS.Surv$Age)
  uni_variate_age = coxph(formula = Surv(Time, Status) ~ Age, data = TS.Surv)
  summary(uni_variate_age)
  
  # CMS
  TS.Surv$CMS = factor(TS.Surv$CMS, levels = c("CMS4", "CMS3", "CMS2", "CMS1"))
  uni_variate_CMS = coxph(formula = Surv(Time, Status) ~ CMS, data = TS.Surv)
  summary(uni_variate_CMS)
  summary =summary(uni_variate_CMS)
  summary$logtest
  
  # Multivariate all
  cox1 = coxph(formula = Surv(Time, Status) ~ ICRscore + MSI + pathologic_stage + Age + CMS, data = TS.Surv)
  summary(cox1)
  
  # Multivariate without cms
  TS.Surv = TS.Surv[-which(is.na(TS.Surv$CMS)),]
  cox2 = coxph(formula = Surv(Time, Status) ~ ICRscore + MSI + pathologic_stage + Age, data = TS.Surv)
  summary(cox2)
  
  anova(cox2, cox1)
}

TS.Surv$GIE_group = factor(TS.Surv$GIE_group, levels = c("less immunoedited","immunoedited"))
TS.Surv$ICR_ordinal = as.character(TS.Surv$ICR_HML)
table(TS.Surv$ICR_ordinal, exclude = NULL)
TS.Surv$ICR_ordinal[which(TS.Surv$ICR_ordinal == "ICR High")] = "3"
TS.Surv$ICR_ordinal[which(TS.Surv$ICR_ordinal == "ICR Medium")] = "2"
TS.Surv$ICR_ordinal[which(TS.Surv$ICR_ordinal == "ICR Low")] = "1"
TS.Surv$ICR_ordinal = as.numeric(TS.Surv$ICR_ordinal)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)

if(cox == "cox_alterative_ordinal_ICR"){
  # Cox regression 
  uni_variate_ICR = coxph(formula = Surv(Time, Status) ~ ICR_ordinal, data = TS.Surv)
  summary(uni_variate_ICR)
  
  # Cox regression MSI
  uni_variate_MSI = coxph(formula = Surv(Time, Status) ~ MSI, data = TS.Surv)
  summary(uni_variate_MSI)
  
  # Univariate stage
  TS.Surv$pathologic_stage = as.numeric(TS.Surv$pathologic_stage)
  uni_variate_stage = coxph(formula = Surv(Time, Status) ~ pathologic_stage, data = TS.Surv)
  summary(uni_variate_stage)
  
  # Age
  TS.Surv$Age = as.numeric(TS.Surv$Age)
  uni_variate_age = coxph(formula = Surv(Time, Status) ~ Age, data = TS.Surv)
  summary(uni_variate_age)
  
  # CMS
  TS.Surv$CMS = factor(TS.Surv$CMS, levels = c("CMS4", "CMS3", "CMS2", "CMS1"))
  uni_variate_CMS = coxph(formula = Surv(Time, Status) ~ CMS, data = TS.Surv)
  summary(uni_variate_CMS)
  summary =summary(uni_variate_CMS)
  summary$logtest
  
  # Multivariate all
  cox1 = coxph(formula = Surv(Time, Status) ~ ICR_ordinal + MSI + pathologic_stage + Age + CMS, data = TS.Surv)
  summary(cox1)
  
  # Multivariate without cms
  TS.Surv = TS.Surv[-which(is.na(TS.Surv$CMS)),]
  cox2 = coxph(formula = Surv(Time, Status) ~ ICR_ordinal + MSI + pathologic_stage + Age, data = TS.Surv)
  summary(cox2)
  
  anova(cox2, cox1)
}

