
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("survival", "plyr", "raster", "texreg", "stringr", "survminer")
ipak(required.packages)

source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

# Set Parameters
Surv.cutoff.years = 20                                        # SET cut-off
Group.of.interest = "Mutational_load_category"  # "Neoantigen_load_category" # "Mutational_load_category"                          # "tumour_anatomic_site" or "ICR_cluster_k3"
CMS = ""
ICR = "ICR Medium"
Survival_outcome = "OS"

# Create folders and log file
dir.create("./Figures/",showWarnings = FALSE)
dir.create("./Figures/WES/006_Mutational_load_KM", showWarnings = FALSE)
dir.create("./Analysis/WES/006_Survival Analysis", showWarnings = FALSE)

# Read in the clinical data file
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

# Divide in groups
frequency_df$Mutational_load_category = NA
frequency_df$Mutational_load_category[which(frequency_df$Nonsilent_mutational_burden_per_Mb < 12)] = "Low"
frequency_df$Mutational_load_category[which(frequency_df$Nonsilent_mutational_burden_per_Mb >= 12)] = "High"

frequency_df$Neoantigen_load_category = NA
frequency_df$Neoantigen_load_category[which(frequency_df$Neoantigen_count < 50)] = "Low"
frequency_df$Neoantigen_load_category[which(frequency_df$Neoantigen_count >= 50)] = "High"


Merged_dataset = clinical_data
Merged_dataset$Mutational_load_category = frequency_df$Mutational_load_category[match(Merged_dataset$Patient_ID,
                                                                                      frequency_df$Patient_ID)]
Merged_dataset$Mutational_load_per_Mb = frequency_df$Nonsilent_mutational_burden_per_Mb[match(Merged_dataset$Patient_ID,
                                                                                              frequency_df$Patient_ID)]
Merged_dataset$Neoantigen_load_category = frequency_df$Neoantigen_load_category[match(Merged_dataset$Patient_ID,
                                                                                      frequency_df$Patient_ID)]
Merged_dataset = Merged_dataset[-which(is.na(Merged_dataset$Mutational_load_category)),]
Merged_dataset$Mutational_load_category = factor(Merged_dataset$Mutational_load_category, levels = c("High", "Low"))
Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                    substring(rownames(table_cluster_assignment), 1, 3))]

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}

if(ICR == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ICR_cluster == ICR),]
}

Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Dead")] = "Event"
Merged_dataset$OS.Status[which(Merged_dataset$OS.Status == "Alive")] = "Event free"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Event")] = "Event"
Merged_dataset$DFS.Status[which(Merged_dataset$DFS.Status == "Disease Free")] = "Event free"

Status = paste0(Survival_outcome, ".Status")
Time = paste0(Survival_outcome, ".Time")

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset[, Status] == "Event free", c(Status, Time, Group.of.interest, "Mutational_load_per_Mb")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "Mutational_load_per_Mb")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset[, Status] == "Event", c(Status, Time, Group.of.interest, "Mutational_load_per_Mb")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "Mutational_load_per_Mb")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Event free"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Event"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

# survival curve
TS.Surv$Mutational_load_category = factor(TS.Surv$Mutational_load_category, levels = c("Low", "High"))
TS.Surv$Neoantigen_load_category = factor(TS.Surv$Neoantigen_load_category, levels = c("Low", "High"))
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv[,Group.of.interest],conf.type = "log-log")

summary = summary(coxph(formula = msurv ~ TS.Surv[,Group.of.interest], data = TS.Surv))
summary
summary$logtest

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("High", "Low"))
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
                         palette = c("#008000", "#9400D3")
)

dir.create("./Figures/WES/006_Mutational_load_KM/006_ggsurv_plots", showWarnings = FALSE)
pdf(paste0("./Figures/WES/006_Mutational_load_KM/006_ggsurv_plots/EDF6b_", Survival_outcome, "_", Group.of.interest, "_",
            ICR, "_Surv_cutoff_years_", Surv.cutoff.years,".pdf"),
    height=4,width=4.8, onefile = FALSE, family = "ArialMT")  # set filename
print(ggsurv_plot)
dev.off()

# Cox regression (continuous)
TS.Surv$Mutational_load_per_Mb = log(TS.Surv$Mutational_load_per_Mb + 1, 10)
uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ Mutational_load_per_Mb, data = TS.Surv)
summary(uni_variate_ICRscore)

# Multivariate
multivariate_ICR_and_stage = coxph(formula = Surv(Time, Status) ~ ICRscore + pathologic_stage, data = TS.Surv)
summary(multivariate_ICR_and_stage)


