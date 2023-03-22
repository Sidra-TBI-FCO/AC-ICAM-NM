
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2","dendextend", "MatrixGenerics"))

load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_patients_microbiome.Rdata")

# subset tumor samples only 
Genus_full_abundance = tcga_bacteria_df

rownames(Genus_full_abundance) = Genus_full_abundance$name
Genus_full_abundance$name = NULL

# create new df
df = data.frame(Genus = rownames(Genus_full_abundance), Number_positive_samples = NA, percentage_positive_samples = NA, maximum_relative_abundance = NA)

# taking maximum value from the row 
Genus_full_abundance = as.matrix(Genus_full_abundance)

df$maximum_relative_abundance = rowMaxs(Genus_full_abundance)

# logical matrix 
genus.new = Genus_full_abundance > 0 

#rowSums(genus.new)

df$Number_positive_samples = rowSums(genus.new)
df$percentage_positive_samples = df$Number_positive_samples/117*100

# skip if per == no
df = df[which(df$percentage_positive_samples >= 10),]
dim(df) # 27 genera remaining


df = df[which(df$maximum_relative_abundance >= 0.01),]

dim(df) # 27 genera remaining

Genus_full_abundance = Genus_full_abundance[df$Genus,]
dim(Genus_full_abundance) # 27 117

tcga_bacteria_df = Genus_full_abundance

save(tcga_bacteria_df, file = "../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_filtered_per_10_RA_0.01.Rdata")
