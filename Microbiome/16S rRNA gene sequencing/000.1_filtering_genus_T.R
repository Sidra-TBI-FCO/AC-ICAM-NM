
## filtering genus based on tumor samples only and create mirror matrix for normal samples validation tumor samples

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2","dendextend", "MatrixGenerics"))

# Set parameters
Type = "Relative"  
Rank = "Genus_full" 
per = 10 # 20 # 10 # 5 # no
abundance = "yes_0.01" # no  yes_0.01 
script = "01.0"

# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))

# subset normal samples (later in script used for mirror)
N_Genus_full_abundance = Genus_full_abundance[,grep("N", colnames(Genus_full_abundance))]

# subset tumor samples only 
Genus_full_abundance = Genus_full_abundance[,grep("T", colnames(Genus_full_abundance))]

# create new df
df = data.frame(Genus = rownames(Genus_full_abundance), Number_positive_samples = NA, percentage_positive_samples = NA, maximum_relative_abundance = NA)

# taking maximum value from the row 
Genus_full_abundance = as.matrix(Genus_full_abundance)

df$maximum_relative_abundance = rowMaxs(Genus_full_abundance)

# logical matrix 
genus.new = Genus_full_abundance > 0 

#rowSums(genus.new)

df$Number_positive_samples = rowSums(genus.new)
df$percentage_positive_samples = df$Number_positive_samples/246*100

# skip if per == no
df = df[which(df$percentage_positive_samples >= per),]
dim(df) # 166 genera remaining

if(abundance == "yes_0.01"){
  df = df[which(df$maximum_relative_abundance >= 0.01),]
}

dim(df) # 138 genera remaining

Genus_full_abundance = Genus_full_abundance[df$Genus,]
dim(Genus_full_abundance) # 138 246

samples_246_T = colnames(Genus_full_abundance)

#dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246")
save(Genus_full_abundance, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/",script,"_Genus_full_abundance_filtered_",
                                         per,"_percent_abundance_",abundance,"_T_246_samples.Rdata"))


################################################################################
# Mirror for 246 Normal samples
N_Genus_full_abundance = N_Genus_full_abundance[df$Genus,]
dim(N_Genus_full_abundance)
# 138 246
Genus_full_abundance = N_Genus_full_abundance

samples_246_N = colnames(Genus_full_abundance)

save(Genus_full_abundance, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/", script, "_Genus_full_abundance_filtered_",
                                         per,"_percent_abundance_",abundance,"_N_246_samples.Rdata"))


################################################################################
# Mirror for 42 Validation tumors
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_all_samples_rank_level_abundancies.Rdata")
Genus_full_abundance = Genus_full_abundance[,-which(colnames(Genus_full_abundance) %in% c(samples_246_T, samples_246_N))]

dim(Genus_full_abundance)
# 562 42 (42 samples)

Genus_full_abundance = Genus_full_abundance[df$Genus,]
dim(Genus_full_abundance)
# 138 42

save(Genus_full_abundance, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/", script, "_Genus_full_abundance_filtered_",
                                         per,"_percent_abundance_",abundance,"_T_42_samples.Rdata"))


Phylum_abundance = Phylum_abundance[,-grep("N", colnames(Phylum_abundance))]
