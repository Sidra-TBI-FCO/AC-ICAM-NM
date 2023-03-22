
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr", "tidyverse")
ipak(required.packages)

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")

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
#matrix = matrix[,c(1:100)]
matrix = matrix[,c(1:200)] # which(colnames(matrix) == "POLE")
  
i = 1
for(i in 1:ncol(matrix)){
  gene = colnames(matrix)[i]
  sub = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  samples = sub$Sample_ID
  matrix[,gene][which(rownames(matrix) %in% samples)] = "MUT"
}

mat = t(matrix)

dir.create("./Analysis/WES/008_matrix_for_Oncoplot", showWarnings = FALSE)
save(mat, file = "./Analysis/WES/008_matrix_for_Oncoplot/matrix_selected_for_oncoplot.Rdata")

