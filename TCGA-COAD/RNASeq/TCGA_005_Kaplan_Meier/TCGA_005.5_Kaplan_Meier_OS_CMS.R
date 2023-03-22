

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr", "survminer")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 20                                        # SET cut-off
Group.of.interest = "CMS"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
exclude_stage_IV = ""
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
Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 12))]
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("CMS1", "CMS2", "CMS3", "CMS4"))

Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$CMS)),]

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

TS.Surv[,Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
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
                         fontsize = 4.5,
                         font.x = 18,
                         font.tickslab = 18,
                         censor.shape = 3,
                         censor.size = 1.5,
                         #pval = TRUE,
                         palette = c("#FF9F21", "#0074AF", "#E97AA8", "#009E74")
)

print(ggsurv_plot)

png(paste0("./Figures/Kaplan Meier Plots/", download.method,"_", exclude_stage_IV ,"_Overall_Survival_CMS_Surv_cutoff_years_", Surv.cutoff.years,".png"),
    res=600,height=3.8,width=4.2,unit="in")  # set filename
print(ggsurv_plot)
dev.off()

# Cox regression
TS.Surv$CMS = factor(TS.Surv$CMS, levels = c("CMS4", "CMS1", "CMS2", "CMS3"))
uni_variate = coxph(formula = Surv(Time, Status) ~ CMS, data = TS.Surv)
sum = summary(uni_variate)
sum
sum$logtest


