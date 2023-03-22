# normal tumor clinical data

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# load data
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
library(dplyr)

tumor = clinical_data
tumor$ID = paste0(tumor$Patient_ID, "T")
tumor <- tumor %>% relocate(ID, .before = gender)

normal = clinical_data
normal$ID = paste0(normal$Patient_ID, "N")
normal <- normal %>% relocate(ID, .before = gender)

clinical_data = rbind(normal, tumor)

clinical_data = clinical_data[order(clinical_data$ID),]

clinical_data$Type = NA
clinical_data$Type[grep("N", clinical_data$ID)] = "Normal"
clinical_data$Type[grep("T", clinical_data$ID)] = "Tumor"

clinical_data <- clinical_data %>% relocate(Type, .before = gender)

save(clinical_data, file = "./Processed_Data/Survival Data/JSREP_NT_clinical_data.Rdata")
