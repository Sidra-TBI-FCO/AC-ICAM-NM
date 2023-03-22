
## Combine data files of TCRSeq kits 1, 2, 3

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Set parameters
data_type = "raw" # "raw" "DNA_input_corrected"

# Load data
TCR_Overview_kit1 = read.csv("./Processed_Data/TCR/TCR SampleOverview/SampleOverview_04-17-2019_10-56-09_PM.csv", stringsAsFactors = FALSE)
TCR_Overview_kit2 = read.csv("./Processed_Data/TCR/TCR SampleOverview/SampleOverview_07-25-2019_6-03-13_AM.csv", stringsAsFactors = FALSE)
TCR_Overview_kit3 = read.csv("./Processed_Data/TCR/TCR SampleOverview/SampleOverview_07-25-2019_5-59-59_AM.csv", stringsAsFactors = FALSE)

# Merge
TCR_Overview_kit2 = TCR_Overview_kit2[-which(TCR_Overview_kit2$sample_name == "114T"),]
TCR_Overview = rbind(TCR_Overview_kit1, TCR_Overview_kit2, TCR_Overview_kit3)
TCR_Overview = TCR_Overview[-which(TCR_Overview$sample_name == "neg"),]

write.csv(TCR_Overview, file = "./Processed_Data/TCR/TCR SampleOverview/SampleOverview_kit_1_2_3.csv")
