
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
selection = "QGPC_genes" # "cancer_driver_501" or "cancer_driver_QGP" or ""
subset = "all_patients"  # "all_patients" or "nonhypermutated" or "hypermutated"

# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load(paste0("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata"))
load("./Processed_Data/External/Gene_collections/501_genes_Michele.Rdata")
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")

finalMafFiltered$Hugo_Symbol = as.character(finalMafFiltered$Hugo_Symbol)
all_patients = unique(finalMafFiltered$Patient_ID)
frequency_df$Patient_ID = as.character(frequency_df$Patient_ID)
nonhypermutated = frequency_df$Patient_ID[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)]
hypermutated = frequency_df$Patient_ID[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)]
table_cluster_assignment$Patient_ID = substring(rownames(table_cluster_assignment), 1, 3)

frequency_df$ICRscore = table_cluster_assignment$ICRscore[match(frequency_df$Patient_ID,
                                                                table_cluster_assignment$Patient_ID)]
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% get(subset)),]
ICR_sub = table_cluster_assignment[which(table_cluster_assignment$Patient_ID %in% get(subset)),]

all_genes = as.character(unique(finalMafFiltered$Hugo_Symbol))
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol %in% get(selection)),]

genes = as.character(unique(finalMafFiltered$Hugo_Symbol))
results = data.frame(Gene = genes, N_mut = NA, 
                     ICR_pval = NA, ICR_FDR = NA, ICR_estimate = NA)
                     #Mutational_load_pval = NA, Mutational_load_estimate = NA)

frequency_df = frequency_df[which(frequency_df$Patient_ID %in% get(subset)),]

frequency_df$histology = clinical_data$Tumor_morphology[match(frequency_df$Patient_ID, clinical_data$Patient_ID)]
#frequency_df = frequency_df[-which(frequency_df$histology == "mucineus adenocarcinoom"),]

#genes = genes[1:10]
i = 1
for (i in 1:length(genes)){
  gene = genes[i]
  sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  df = frequency_df[, c("Patient_ID", "ICRscore", "Non_silent_Mutation_frequency")]
  df$Mutation = 0
  patients = unique(sub$Patient_ID)
  df$Mutation[which(df$Patient_ID %in% patients)] = 1
  
  #fit = glm(Mutation ~ ICRscore + Non_silent_Mutation_frequency, data = df, family = "binomial") with covariate
  fit = glm(Mutation ~ ICRscore, data = df, family = "binomial")
  table = summary(fit)
  test_table = data.frame(table$coefficients)
  results$N_mut[which(results$Gene == gene)] = sum(df$Mutation)
  results$ICR_pval[which(results$Gene == gene)] = test_table["ICRscore","Pr...z.."]
  results$ICR_estimate[which(results$Gene == gene)] = test_table["ICRscore","Estimate"]
  #results$Mutational_load_pval[which(results$Gene == gene)] = test_table["Non_silent_Mutation_frequency","Pr...z.."]
  #results$Mutational_load_estimate[which(results$Gene == gene)] = test_table["Non_silent_Mutation_frequency","Estimate"]
}
results = results[order(results$ICR_pval),]
results$ICR_FDR = p.adjust(results$ICR_pval, method = "BH", n = nrow(results))

dir.create("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR", showWarnings = FALSE)
write.csv(results, file = paste0("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR/", subset,
                                 "_results_linear_regression_model_Somatic_mutation_by_ICRscore.csv"),
          row.names = FALSE)

results = results[-which(results$N_mut <10),]
results$ICR_FDR = p.adjust(results$ICR_pval, method = "BH", n = nrow(results))

write.csv(results, file = paste0("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR/Min_10_Mut_", subset,
                                 "_results_linear_regression_model_Somatic_mutation_by_ICRscore.csv"),
          row.names = FALSE)
