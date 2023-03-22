
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

dir.create("./Analysis/WES", showWarnings = FALSE)

# Load data
load("./Processed_Data/WES/MAF/finalMaf_clean_primary_tumor_samples.Rdata")
load("./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_clean_samples.Rdata")

# All mutation rate (silent and non-silent)
table(finalMaf$Variant_Classification)
table(finalMaf$Tissue)
table(finalMafFiltered$Tissue)

frequency_df = data.frame(table(finalMaf$Patient_ID))
colnames(frequency_df) = c("Patient_ID", "Mutation_frequency")

frequency_df$Total_mutational_burden_per_Mb = frequency_df$Mutation_frequency / 40

dir.create("./Analysis/WES/002_Mutation_frequency", showWarnings = FALSE)
save(frequency_df, file = paste0("./Analysis/WES/002_Mutation_frequency/Total_mutation_frequency.Rdata"))


# 
table(finalMafFiltered$Variant_Classification)

frequency_df = data.frame(table(finalMafFiltered$Patient_ID))
colnames(frequency_df) = c("Patient_ID", "Non_silent_Mutation_frequency")

frequency_df$Nonsilent_mutational_burden_per_Mb = frequency_df$Non_silent_Mutation_frequency / 40

dir.create("./Analysis/WES/002_Mutation_frequency", showWarnings = FALSE)
save(frequency_df, file = paste0("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata"))

