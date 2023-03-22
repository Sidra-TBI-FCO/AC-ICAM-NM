
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak("stringr")

# Load data
Training = read.csv("./Processed_Data/Microbiome/External/9_August/Version2_ACICAM_Risk_Scores.csv", stringsAsFactors = FALSE)
Validation = read.csv("./Processed_Data/Microbiome/External/9_August/Version2_ACICAM_Subset_42_Risk_Scores.csv", stringsAsFactors = FALSE)

Training$Patient_ID = str_pad(pad = "0", width =3, Training$samples)
Validation$Patient_ID = str_pad(pad = "0", width =3, Validation$samples)

save(Training, Validation, file = "./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
