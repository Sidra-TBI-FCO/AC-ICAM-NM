

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 10                                        # SET cut-off
Group.of.interest = "ICR_HML"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_stage_IV = ""
CMS = ""
download.method = "Biolinks"

# Create folders and log file
dir.create("./Figures",showWarnings = FALSE)
dir.create("./Figures/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Analysis", showWarnings = FALSE)
dir.create("./Analysis/Survival Analysis", showWarnings = FALSE)

# Read in the clinical data file
clinical_data = read.csv("./Processed_Data/External_Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                         stringsAsFactors = FALSE)
colnames(clinical_data)[which(colnames(clinical_data) == "bcr_patient_barcode")] = "Patient_ID"

# Add ICR as a variable and assign ICR cluster according to table cluster assignment
if(download.method == "Biolinks"){
  load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
  load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
  load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
  table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(filtered.norm.RNAseqData)),]
}
if(download.method == "TCGA_Assembler"){
  load("./Analysis/ICR Consensus Clustering/TCGA_Assembler_COAD_ICR_cluster_assignment_k2-6.Rdata")
  load("./Processed_Data/RNASeq/TCGA_Assembler_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
  load("./Analysis/016_CMS_Classification/TCGA_Assembler_Rfcms.Rdata")
}

colnames(table_cluster_assignment)[which(colnames(table_cluster_assignment) == "HML_cluster")] = "ICR_HML"

table_cluster_assignment$Patient_ID = substring(row.names(table_cluster_assignment), 1, 12)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("ICR High", "ICR Medium", "ICR Low"))
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

# Exclude adjuvant treated
# not applicable for TCGA

# Exclude Stage IV
if(exclude_stage_IV == "exclude_stage_IV" ){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ajcc_pathologic_tumor_stage == "IV"),]
}

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}


# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$vital_status == "Alive", c("vital_status", "last_contact_days_to", Group.of.interest, "ICRscore")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "ICRscore")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$vital_status == "Dead", c("vital_status", "death_days_to", Group.of.interest, "ICRscore")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "ICRscore")
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

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
#Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
#Cal.Surv[,Group.of.interest] = relevel(Cal.Surv[,Group.of.interest], "ICR3")
Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("ICR High", "ICR Medium", "ICR Low"))
summary(mHR)
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)


# plots
png(paste0("./Figures/Kaplan Meier Plots/", download.method, "_", exclude_stage_IV, "_", CMS,"_Overall_Survival_", Group.of.interest,"_Surv_cutoff_years_",Surv.cutoff.years,".png"),res=600,height=6,width=8,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     palette = c("red", "green", "blue"),
     ylabs = "OS probability",
     PLOT_P = signif(p[2],3),
     PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3),
     PLOT_CI1 = CI[2,1],
     PLOT_CI2 = CI[2,2])
dev.off()

dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,Group.of.interest]),
     ystrataname = "Legend",
     main= paste0("Survival curve across ", Group.of.interest),
     xlabs = "Time in months",
     palette = c("red", "green", "blue"),
     ylabs = "OS probability",
     #PLOT_P = signif(p[2],3),
     #PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3),
     #PLOT_CI1 = CI[2,1],
     #PLOT_CI2 = CI[2,2],
     legend = FALSE)

# Cox regression (continuous)
uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
summary(uni_variate_ICRscore)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)

