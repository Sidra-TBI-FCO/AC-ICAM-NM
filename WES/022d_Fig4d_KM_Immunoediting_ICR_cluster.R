

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))

ipak(c("forestplot", "survminer", "survival"))

# Set parameters
exclude_medium = "exclude_medium" # "include_medium_cat" "include_medium" or "exclude_medium"
Group.of.interest = "Immunoedited_ICR_cluster"
Surv.cutoff.years = 20
x = 2
Stages = "all stages" #"all stages"  #c(1,2,3) #"all stages" #"all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name =  "All" #"StageI-II-III" #"Stage III" # "All" #"Stage III&IV" #"All"
CMS = "" # "" for all, "CMS1" or "CMS2"

# Load data
load(paste0("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset_",
            exclude_medium, ".Rdata"))

if(exclude_medium == "include_medium_cat"){
  #palette = c("#FF3805", "#FF9F05", "#BEBEBE", "#BBBBBB", "#039EFC", "#5203FC")
  palette = c("#FF3805", "#FF3806", "#4EAD5B", "#00B051", "#0030FF", "#0030FE")
  levels1 = c("less immunoedited ICR Low", "immunoedited ICR High", "less immunoedited ICR High",
             "immunoedited ICR Medium", "less immunoedited ICR Medium",
             "immunoedited ICR Low")
  levels2 = c("immunoedited ICR High", "less immunoedited ICR High",
              "immunoedited ICR Medium", "less immunoedited ICR Medium",
              "immunoedited ICR Low", "less immunoedited ICR Low")
  linetypes = c("solid","dashed","solid","dashed", "solid", "dashed")
}else{
  palette = c("#FF3805", "#FF9F05", "#039EFC", "#5203FC")
  levels1 = c("less immunoedited ICR Low", "immunoedited ICR High", "less immunoedited ICR High",
             "immunoedited ICR Low")
  levels2 = c("immunoedited ICR High", "less immunoedited ICR High",
              "immunoedited ICR Low", "less immunoedited ICR Low")
  linetypes = c("solid","solid","solid","solid", "solid", "solid")
}

# Overview subgroups
table(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)
Merged_dataset$Immunoedited_ICR_cluster = paste(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)

table(Merged_dataset$Immunoedited_ICR_cluster)
table(Merged_dataset$ajcc_pathologic_tumor_stage)
if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}
table(Merged_dataset$ajcc_pathologic_tumor_stage)

if(CMS == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$CMS == CMS),]
}
table(Merged_dataset$CMS)

### KM
# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Merged_dataset[Merged_dataset$OS.Status == "Alive", c("OS.Status", "OS.Time", Group.of.interest,
                                                                 "Immunoediting_score", "Nonsilent_mutational_burden_per_Mb")]
colnames(TS.Alive) = c("Status","Time", Group.of.interest, 
                       "Immunoediting_score", "Nonsilent_mutational_burden_per_Mb")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Merged_dataset[Merged_dataset$OS.Status == "Dead", c("OS.Status", "OS.Time", Group.of.interest,
                                                               "Immunoediting_score", "Nonsilent_mutational_burden_per_Mb")]
colnames(TS.Dead) = c("Status","Time", Group.of.interest,
                      "Immunoediting_score", "Nonsilent_mutational_burden_per_Mb")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "Alive"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "Dead"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

TS.Surv[, Group.of.interest] = factor(TS.Surv[,Group.of.interest], levels =  levels2)
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
Cal.Surv[,Group.of.interest] = factor(Cal.Surv[,Group.of.interest], levels = 
                                        levels1)


mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv, subset = Cal.Surv[, Group.of.interest] %in% c("immunoedited ICR High", "less immunoedited ICR Low"))

mHR = coxph(formula = msurv ~ Cal.Surv[,Group.of.interest],data = Cal.Surv)
summary = summary(mHR)
summary$logtest
summary

Cal.Surv$Immunoedited_ICR_cluster_ordinal = Cal.Surv$Immunoedited_ICR_cluster
#levels(Cal.Surv$Immunoedited_ICR_cluster_ordinal) = c("1", "3", "2", "2")
levels(Cal.Surv$Immunoedited_ICR_cluster_ordinal) = c("1", "2", "2", "2")
Cal.Surv$Immunoedited_ICR_cluster_ordinal = as.character(Cal.Surv$Immunoedited_ICR_cluster_ordinal)
Cal.Surv$Immunoedited_ICR_cluster_ordinal = as.numeric(Cal.Surv$Immunoedited_ICR_cluster_ordinal)

ordinal_coxph = coxph(formula = msurv ~ Cal.Surv$Immunoedited_ICR_cluster_ordinal, data = Cal.Surv)
summary = summary(ordinal_coxph)
summary$logtest
summary

# the copxph of scale(ICR score) + scale(1/IE)
Cal.Surv$continuous_ICR_IE_score 

# coxph using tmb and immunoediting as continuous
Cal.Surv$Mutation_rate = log10(Cal.Surv$Nonsilent_mutational_burden_per_Mb + 1)
multivariate = coxph(formula = Surv(Time, Status) ~ Immunoediting_score + Mutation_rate, data = Cal.Surv)
summary(multivariate)

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
                         linetype = linetypes,
                         #pval = TRUE,
                         palette = palette
)

dir.create("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/ggsurvplots", showWarnings = FALSE)
pdf(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/ggsurvplots/Fig4d_ggsurv_Dec_2020_OS_by_ICR_cluster_and_immunoediting_", Stages, "_",
           Group.of.interest, exclude_medium, CMS, ".pdf"),
    #res = 600, height = 6, width = 5.4, units = "in")
    #res = 600, height = 4.2, width = 5.4, units = "in")
    height=6,width=5.4, onefile = FALSE, family = "ArialMT")
print(ggsurv_plot)
dev.off()
