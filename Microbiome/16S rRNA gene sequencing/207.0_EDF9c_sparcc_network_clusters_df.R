
# clusters df v1:

# Sum up the RA: Include only G1, 2, 3, and 5. For G4, exclude Prevotella 2, the metrics of the other Gs (1,2,3, and 5) will remain the same as there are not negatively correlated taxa

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx", "stringr"))

# Set parameters
per = 10 # 10 # 5 no
abundance = "yes_0.01" # no  yes_0.01 
clusters.v1 = c("Group1", "Group2", "Group3", "Group4", "Group5")
samples = "tumors" # tumors # normal
matrix = "RA" # RA # OTU

# load data
network = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/network/sparcc_network_analysis_genera_v3.xlsx")
if (matrix == "RA"){
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
}

if (matrix == "OTU"){
  load("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/matrices/OTU_sparcc_matrix_genera_names_fixed.Rdata")
  Genus_full_abundance = sparcc.fil
}

network$clutser = str_trim(network$clutser, "both")

# exclude groups 6 and 7
network = network[which(network$clutser %in% clusters.v1),]
table(network$clutser)

# exclude prevotella 2  
network = network[-grep("Prevotella 2", network$genera),]
table(network$clutser)

df = data.frame(paitent_ID = colnames(Genus_full_abundance), Group1 = NA, Group2 = NA, Group3 = NA, Group4 = NA, Group5 = NA)

Genus_full_abundance = as.data.frame(Genus_full_abundance)
Genus_full_abundance$genus = gsub(".*\\D_5__", "", rownames(Genus_full_abundance))

network$genera[which(network$genera %in% Genus_full_abundance$genus)]

i = 1
for (i in 1:length(clusters.v1)) {
  cluster = clusters.v1[i]
  genus = network[which(network$clutser == cluster),]
  relative.abundance = Genus_full_abundance[which(Genus_full_abundance$genus %in% genus$genera),]
  relative.abundance$genus = NULL
  relative.abundance = as.matrix(relative.abundance)
  df[,cluster] = colSums(relative.abundance)
}

save(df, file = paste0("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/network/network_clusters_df_without_Prevotella2_",matrix,".Rdata"))


