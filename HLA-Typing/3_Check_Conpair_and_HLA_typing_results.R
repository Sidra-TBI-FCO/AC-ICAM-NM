

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr")
ipak(required.packages)

# Load data
discon_details_df = read.csv("./Analysis/HLA Typing/Optitype_concordance/disconcordancy_details.csv", stringsAsFactors = FALSE)
conpair_result = read.csv("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/ConPair/Conpair_result.csv", stringsAsFactors = FALSE)

# Analysis
conpair_result$Patient_ID = gsub("N.*", "", conpair_result$Sample_combination)

Patients_1_mismatch_TWNW = discon_details_df$patient_ID[which(discon_details_df$TWNW.mismatch == 1)]
Patients_high_conpair = conpair_result$Patient_ID[which(conpair_result$Concordance >0.98)]

length(intersect(Patients_1_mismatch_TWNW, Patients_high_conpair))
