
### Correlation genes and clonality

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
test = "pearson"

# load data
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")

# subset
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
dim(RNASeq.QN.LOG2)
RNA_expression = RNASeq.QN.LOG2[, which(substring(colnames(RNASeq.QN.LOG2), 1, 4) %in% TCR_Overview$sample_name)]
dim(RNA_expression)

# remove genes with all values 0
RNA_expression = RNA_expression[-which(rowSums(RNA_expression) == 0),]

colnames(RNA_expression) = substring(colnames(RNA_expression), 1, 4)

df = data.frame(t(RNA_expression))
df$productive_clonality = TCR_Overview$productive_clonality[match(rownames(df), TCR_Overview$sample_name)]

results = data.frame(Gene = colnames(df)[1:(ncol(df)-1)],
                     Corr_coefficient = NA,
                     p_value = NA,
                     CI_lower = NA,
                     CI_upper = NA)
results$Gene = as.character(results$Gene)

i=1  
for (i in 1:nrow(RNA_expression)){
  gene = results$Gene[i]
  cor_test = cor.test(df[, gene], df$productive_clonality, method = test, conf.level = 0.95)
  results$Corr_coefficient[which(results$Gene == gene)] = cor_test$estimate
  results$p_value[which(results$Gene == gene)] = cor_test$p.value
  results$CI_lower[which(results$Gene == gene)] = cor_test$conf.int[1]
  results$CI_upper[which(results$Gene == gene)] =  cor_test$conf.int[2]
}

dir.create("./Analysis/TCR/009_Correlation_gene_expression_clonality", showWarnings = FALSE)
write.csv(df, file = paste0("./Analysis/TCR/009_Correlation_gene_expression_clonality/114_all_vars_df.csv"))
write.csv(results, file = paste0("./Analysis/TCR/009_Correlation_gene_expression_clonality/114_",
                                 test, "_correlation_genes_productive_clonality.csv"),
          row.names = FALSE)
