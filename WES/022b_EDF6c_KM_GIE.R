
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(required.packages = c("survival", "plyr", "survminer"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

# Set parameters
cutoff = 1
subset = ""  #"nonhypermutated" "hypermutated"
Group.of.interest = "GIE"  #"GIE"
x = 1
Surv.cutoff.years = 20
Stages = "all stages" #"all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name = "All" #"Stage I&II" #"All" # "Stage III&IV"
ICR_filter = "" # "" or "ICR Medium"

# Load data
load(paste0("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata"))
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Survival analysis
IES_df$GIE = gsub(" ICR*.*", "", IES_df$Immunoedited_ICR_cluster)
Merged_dataset = clinical_data[which(clinical_data$Patient_ID %in% IES_df$Patient_ID),]
Merged_dataset = merge(Merged_dataset, IES_df, by = "Patient_ID")

table(Merged_dataset$ajcc_pathologic_tumor_stage)
if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}
table(Merged_dataset$ajcc_pathologic_tumor_stage)

Merged_dataset$ajcc_pathologic_tumor_stage = as.numeric(Merged_dataset$ajcc_pathologic_tumor_stage)
Merged_dataset$age_at_initial_pathologic_diagnosis = as.numeric(Merged_dataset$age_at_initial_pathologic_diagnosis)
Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                    substring(rownames(Merged_dataset), 1, 3))]

if(ICR_filter == ""){}else{
  Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                      substring(rownames(table_cluster_assignment), 1, 3))]
  Merged_dataset = Merged_dataset[which(Merged_dataset$ICR_cluster == ICR_filter),]
  
}

### KM
# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$OS.Status == "Alive", c("OS.Status", "OS.Time", Group.of.interest,
                                                                 "ajcc_pathologic_tumor_stage", 
                                                                 "age_at_initial_pathologic_diagnosis")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, "Stage", "Age")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$OS.Status == "Dead", c("OS.Status", "OS.Time", Group.of.interest,
                                                              "ajcc_pathologic_tumor_stage",
                                                               "age_at_initial_pathologic_diagnosis")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest, "Stage", "Age")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

TS.Surv[,Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels = c("less immunoedited", "immunoedited"))

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
if(Group.of.interest == "Immunoedited_extremes"){
  Cal.Surv[,Group.of.interest] = factor(Cal.Surv[,Group.of.interest], levels = c("immunoedited", "medium immunoediting", "less immunoedited"))
}

mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("immunoedited", "less immunoedited"))
summary = summary(mHR)
summary$logtest
summary

if(x == 1){
  palette = c("darkblue", "orange")
}
if(x == 2){
  palette = c("#FF5800", "#FFCC00","#8085E9")
}

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
                         pval = TRUE,
                         palette = palette)

dir.create("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/ggsurvplots", showWarnings = FALSE)
pdf(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/ggsurvplots/EDF6c_", cutoff, "_ggsurv_", subset, "_by_immunoediting_status_", 
           Group.of.interest, "_", Stage_name, "_",ICR_filter,".pdf"),
    height=4,width=4.8, onefile = FALSE, family = "ArialMT")
print(ggsurv_plot)
dev.off()

summary
summary$logtest
