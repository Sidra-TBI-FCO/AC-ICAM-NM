
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 20                                        # SET cut-off
Group.of.interest = "CMS"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
ICR = ""
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Create folders and log file
dir.create("./Figures/Trimmed_p",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Logfiles/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p", showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p/Survival Analysis", showWarnings = FALSE)

# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Add ICR as a variable and assign ICR cluster according to table cluster assignment
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
table_cluster_assignment$Patient_ID = row.names(table_cluster_assignment)
table_cluster_assignment$Patient_ID = gsub("T_P", "", table_cluster_assignment$Patient_ID)

Merged_dataset = merge(clinical_data, table_cluster_assignment, by = "Patient_ID")

Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]
Merged_dataset[, Group.of.interest] = factor(Merged_dataset[, Group.of.interest], levels = c("CMS4", "CMS3", "CMS2", "CMS1"))

Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$CMS)),]

if(exclude == ""){
}else{
  Merged_dataset = Merged_dataset[-which(Merged_dataset$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
}

if(ICR > 1){
  Merged_dataset = Merged_dataset[which(Merged_dataset$ICR_HML == ICR),]
}


# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$DFS.Status == "Disease Free", c("DFS.Status", "DFS.Time", Group.of.interest, "ICRscore")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "ICRscore")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$DFS.Status == "Event", c("DFS.Status", "DFS.Time", Group.of.interest, "ICRscore")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "ICRscore")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Disease Free"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Event"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                              # remove patients with less then 1 day follow up time

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
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("CMS1", "CMS2", "CMS3", "CMS4"))
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

summary(mHR)
logrank_test(Surv(Time, Status) ~ get(Group.of.interest), data = TS.Surv[which(TS.Surv[, Group.of.interest] %in% c("CMS1", "CMS4")),])

# plot
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
                         palette = c("#FF9F21", "#0074AF", "#E97AA8", "#009E74")
)

dir.create("./Figures/Trimmed_p/005_ggsurv_plots", showWarnings = FALSE)
pdf(paste0("./Figures/Trimmed_p/005_ggsurv_plots/Fig1e_", ICR, "_", str_c(exclude, collapse = "_"),
           "_PFS_", Group.of.interest,"_Surv_cutoff_years_",Surv.cutoff.years,".pdf"),
    height=4,width=4.8, onefile = FALSE, family = "ArialMT")  # set filename
print(ggsurv_plot)
dev.off()
