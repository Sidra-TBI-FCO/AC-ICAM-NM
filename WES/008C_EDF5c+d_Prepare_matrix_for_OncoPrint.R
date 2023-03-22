

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr", "tidyverse")
ipak(required.packages)

# Set parameters
version = "glm_2"  # "glm_1" for linear model genes that inversely correlate with ICR, "glm_2" for linear model genes that positively correlate with ICR, glm_3 overall

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")
DDR = read.csv("./Processed_Data/External/Gene_collections/IND99_GENES_ANNOTATION_JR_V5_DB.csv", stringsAsFactors = FALSE)
all_patients = read.csv(paste0("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR/all_patients_results_linear_regression_model_Somatic_mutation_by_ICRscore.csv"), stringsAsFactors = FALSE)
hypermutated = read.csv(paste0("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR/hypermutated_results_linear_regression_model_Somatic_mutation_by_ICRscore.csv"), stringsAsFactors = FALSE)
nonhypermutated = read.csv(paste0("./Analysis/WES/025_Linear_Regression_Model_Somatic_Mutation_and_ICR/nonhypermutated_results_linear_regression_model_Somatic_mutation_by_ICRscore.csv"), stringsAsFactors = FALSE)


# Filter minimum number of mutations
if(version %in% c("glm_1")){
  all_patients = all_patients[which(all_patients$ICR_FDR < 0.1 & all_patients$N_mut > 4 & all_patients$ICR_estimate < 0),]
  nonhypermutated = nonhypermutated[which(nonhypermutated$ICR_pval < 0.1 & nonhypermutated$N_mut > 4 & nonhypermutated$ICR_estimate < 0),]
  hypermutated = hypermutated[which(hypermutated$ICR_pval < 0.05 & hypermutated$N_mut > 4 & hypermutated$ICR_estimate < 0),]
}

if(version %in% c("glm_2")){
  all_patients = all_patients[which(all_patients$ICR_FDR < 0.05 & all_patients$N_mut > 4 & all_patients$ICR_estimate > 0),]
  nonhypermutated = nonhypermutated[which(nonhypermutated$ICR_pval < 0.05 & nonhypermutated$N_mut > 4 & nonhypermutated$ICR_estimate > 0),]
  hypermutated = hypermutated[which(hypermutated$ICR_pval < 0.05 & hypermutated$N_mut > 4 & hypermutated$ICR_estimate > 0),]
}

if(version %in% c("glm_3")){
  all_patients = all_patients[which(all_patients$ICR_FDR < 0.05 & all_patients$N_mut > 4),]
}

intersect(nonhypermutated$Gene, hypermutated$Gene)

if(version == "glm_1"){
  nonhypermutated = nonhypermutated[c(1, 3, 2),]
}


if(version %in% c("glm_1", "glm_2")){
  genes = c(nonhypermutated$Gene, hypermutated$Gene)
}

if(version == "glm_3"){
  genes = all_patients$Gene
}

# Gene count ranking
finalMafFiltered$Tumor_Sample_Barcode = as.character(finalMafFiltered$Tumor_Sample_Barcode)
finalMafFiltered$Hugo_Symbol = as.character(finalMafFiltered$Hugo_Symbol)
MAF_gene_Count= finalMafFiltered %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count))

matrix = matrix(nrow = length(unique(finalMafFiltered$Sample_ID)), ncol = nrow(MAF_gene_Count))
colnames(matrix) = MAF_gene_Count$Hugo_Symbol
rownames(matrix) = unique(finalMafFiltered$Sample_ID)

matrix[which(is.na(matrix))] = ""
matrix = matrix[,genes]
#matrix = matrix[,c(148, 1:100)] # 132 = POLE

i = 1
for(i in 1:ncol(matrix)){
  gene = colnames(matrix)[i]
  sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  samples = sub$Sample_ID
  matrix[,gene][which(rownames(matrix) %in% samples)] = "MUT"
}

mat = t(matrix)

dir.create("./Analysis/WES/008_matrix_for_Oncoplot", showWarnings = FALSE)
save(mat, file = paste0("./Analysis/WES/008_matrix_for_Oncoplot/008", 
                        version, "_matrix_for_oncoplot.Rdata"))

