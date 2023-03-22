
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("reshape2", "ggplot2", "ggpubr"))

# Load data
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
TCGA_SScms = SScms
load("../NGS_data/Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
AC_ICAM_SScms = SScms

# 
TCGA_SScms$Cohort = "TCGA-COAD"
AC_ICAM_SScms$Cohort = "AC_ICAM"
TCGA_SScms$SSP.predictedCMS[which(is.na(TCGA_SScms$SSP.predictedCMS))] = "Mixed"
AC_ICAM_SScms$SSP.predictedCMS[which(is.na(AC_ICAM_SScms$SSP.predictedCMS))] = "Mixed"

df = data.frame(prop.table(table(AC_ICAM_SScms$SSP.predictedCMS)))
colnames(df) = c("CMS", "AC_ICAM")
df$TCGA_COAD = NA

df_TCGA = data.frame(prop.table(table(TCGA_SScms$SSP.predictedCMS)))
df$TCGA_COAD = df_TCGA$Freq[match(df$CMS, df_TCGA$Var1)]

df_melt = melt(df, id.vars = "CMS")
df_melt$CMS = factor(df_melt$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "Mixed"))

png("../NGS_data/Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/SSP_CMS_assignment_by_cohort.png",
    res = 600, units = "in", width = 3, height = 3)
ggplot(df_melt, aes(x = variable, y = value, fill = CMS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                               "Mixed" = "lightgrey")) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  ylab("Proportion") +
  xlab("")
dev.off()

pdf("../NGS_data/Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/SSP_CMS_assignment_by_cohort.pdf",
     width = 3, height = 3)
ggplot(df_melt, aes(x = variable, y = value, fill = CMS)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                               "Mixed" = "lightgrey")) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  ylab("Proportion") +
  xlab("")
dev.off()

df_combined = rbind(AC_ICAM_SScms, TCGA_SScms)

table(df_combined$Cohort, df_combined$SSP.predictedCMS)

df_combined$SSP.predictedCMS[-which(df_combined$SSP.predictedCMS == "CMS2")] = "Other"

chisq.test(table(df_combined$Cohort, df_combined$SSP.predictedCMS))
