
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "stringr"))

# Load data
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
model.genera = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv")

model = model.genera$covariates

# subset RA matrix to selected genera
Genus_full_abundance = Genus_full_abundance[which(rownames(Genus_full_abundance) %in% model),]
Genus_full_abundance = t(Genus_full_abundance)
Genus_full_abundance = as.data.frame(Genus_full_abundance)

# add risk score
Training$Patient_ID = paste0(Training$samples, "T")
Training$Patient_ID = str_pad(Training$Patient_ID, 4, pad = 0)

Genus_full_abundance$risk.score = Training$prediction[match(rownames(Genus_full_abundance), Training$Patient_ID)]
Genus_full_abundance = t(Genus_full_abundance)

# Analysis
abundance_T = Genus_full_abundance
abundance_T = t(abundance_T)

# fix rownames 
rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 4)

# Select right rows (patients)
immune_sig_df = immune_sig_df[which(rownames(immune_sig_df) %in% rownames(abundance_T)),]

# Analysis
df_plot = data.frame(Patient_ID = rownames(abundance_T),  Enrichment = NA, Abundance = NA)

results_all = data.frame(Name = NA, signature = NA, p_val = NA, rho = NA)

rank_abundancies = Genus_full_abundance


j = 8 # 22 (= TREM1 data) # 78 (= attractor G) # 8 (PDL1 data) # 63 (GRANS PCA)
for(j in 1:ncol(immune_sig_df)){
  sig = colnames(immune_sig_df)[j]
  df_plot$Enrichment = immune_sig_df[, sig][match(df_plot$Patient_ID, rownames(immune_sig_df))]
  
  results = data.frame(Name = rownames(rank_abundancies), signature = NA, p_val = NA, rho = NA)
  
  i= 371 # 371 (= Ruminococcus 1) # 374 (= subdo) 414 (=Selenomonas) 277 (Fusicatenibacter)
  for (i in 1:nrow(rank_abundancies)){
    micr = rownames(rank_abundancies)[i]
    df_plot$Abundance = abundance_T[, micr][match(df_plot$Patient_ID, rownames(abundance_T))]
    
    df_calc = df_plot
    df_plot_edit = df_plot
    df_calc$Enrichment = as.numeric(df_calc$Enrichment)
    cor = cor.test(df_calc$Enrichment, df_calc$Abundance, method = "pearson") # pearson # spearman 
    results$p_val[which(results$Name == micr)] = cor$p.value
    results$rho[which(results$Name == micr)] = cor$estimate
    
  }
  
  results$signature = sig
  results_all = rbind(results_all, results)
}

results_all = results_all[!is.na(results_all$p_val),]
results_all = results_all[order(results_all$p_val),]

