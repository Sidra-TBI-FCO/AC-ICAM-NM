
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages =  c("stringr", "ggplot2")
ipak(required.packages)

MANTIS = read.csv("./Processed_Data/WES/MANTIS/allsampleMANTIS.csv", stringsAsFactors = FALSE)

MANTIS$Patient_ID = gsub("\\..*", "", MANTIS$BAM.files)
MANTIS$Patient_ID = gsub("N", "", MANTIS$Patient_ID)
MANTIS$Patient_ID = str_pad(MANTIS$Patient_ID , 3, pad = "0") # add leading zero
MANTIS$Sample_ID = gsub("*.sorted.bam_file", "", MANTIS$BAM.files)
MANTIS$Sample_ID = gsub(".*\\_", "", MANTIS$Sample_ID)
MANTIS$Tissue = substring(MANTIS$Sample_ID, nchar(MANTIS$Sample_ID), nchar(MANTIS$Sample_ID))
MANTIS$Tissue[which(MANTIS$Tissue == "M")] = "LM"
MANTIS$Tissue[which(MANTIS$Tissue == "1")] = "LM1"
MANTIS$Tissue[which(MANTIS$Tissue == "2")] = "LM2"
MANTIS$Sample_ID = paste(MANTIS$Patient_ID, MANTIS$Tissue, sep = "")

MANTIS$MSI = NA
MANTIS$MSI[which(MANTIS$MANTIS.score > 0.4)] = "MSI-H"
MANTIS$MSI[which(MANTIS$MANTIS.score <= 0.4)] = "MSS"

table(MANTIS$MSI)

save(MANTIS, file = "./Processed_Data/WES/MANTIS/MANTIS.Rdata")
