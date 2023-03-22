

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr")
ipak(required.packages)

# Load data
HLA_type = read.csv("./Processed_Data/HLA-Typing/WES_HLA_Optitype/Final_HLA_Optitype_results_WES.csv",
                    stringsAsFactors = FALSE)
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, pad = 0, 3)
excluded_df = excluded_df[which(excluded_df$Reason.excluded == "Conpair_lower_90_percent"),]

# Prepare data
HLA_type$Sample_ID = gsub("WES_COAD_LUMC_SIDRA_", "", HLA_type$sample)
HLA_type$Sample_ID = str_sub(HLA_type$Sample_ID, end =-11)
HLA_type$Sample_ID = str_pad(HLA_type$Sample_ID, pad = 0, 4)

HLA_type$Patient_ID = substring(HLA_type$Sample_ID, 1, 3)
HLA_type = HLA_type[-which(HLA_type$Patient_ID %in% excluded_df$Patient_ID),]
HLA_type$Tissue_type = substring(HLA_type$Sample_ID, 4, nchar(HLA_type$Sample_ID))

HLA_type = HLA_type[order(HLA_type$Patient_ID),]

# Check
HLA_type_T = HLA_type[which(HLA_type$Tissue_type == "T"),]
patients = unique(HLA_type_T$Patient_ID)

HLA_type_N = HLA_type[which(HLA_type$Tissue_type == "N"),]
test = HLA_type_N[which(HLA_type_N$Patient_ID %in% patients),]
patients_N = unique(test$Patient_ID)

write.csv(HLA_type, file = "./Processed_Data/HLA-Typing/WES_HLA_Optitype/v2_Filtered_Sidra_LUMC_Cohort_Optitype.csv", row.names = FALSE)
save(HLA_type_T, HLA_type_N, file = "./Processed_Data/HLA-Typing/WES_HLA_Optitype/v2_Filtered_Sidra_LUMC_Cohort_Optitype.Rdata")

