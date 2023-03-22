
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier/HR_table_OS_246.Rdata")
HR_table_246 = results_df

load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier/HR_table_OS_42.Rdata")
HR_table_42 = results_df

load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier/HR_table_OS_TCGA+42_Validation.Rdata")
HR_table_42_TCGA_COAD = results_df

load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/101_Kaplan_Meier/HR_table_OS_TCGA.Rdata")
HR_table_TCGA_COAD = results_df

# Combined all
HR_table_all = rbind(HR_table_246, HR_table_42)
HR_table_all = rbind(HR_table_all, HR_table_TCGA_COAD)
HR_table_all = rbind(HR_table_all, HR_table_42_TCGA_COAD)

dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/105_Collect_Data_for_forest_plot", showWarnings = FALSE)

save(HR_table_all, file = "./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/105_Collect_Data_for_forest_plot/105_HR_table_AC_ICAM_ICAM42_TCGA_COAD.Rdata")
