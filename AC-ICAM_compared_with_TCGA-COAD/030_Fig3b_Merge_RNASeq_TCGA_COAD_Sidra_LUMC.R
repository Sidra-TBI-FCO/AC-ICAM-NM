
# Merge RNASeq TCGA-COAD and Sidra-LUMC and normalize together (EDASeq, QN, log transform)

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("TCGAbiolinks", "EDASeq", "preprocessCore")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters

# Load data
# TCGA-COAD Biolinks
load("~/Dropbox (SMCP)/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica/Processed_Data/RNASeq/001_Exp_COAD_Processed_Data_by_GeneSymbol.rda")
load("~/Dropbox (Personal)/Bioinformatics tools/GeneInfo/geneInfo.July2017.RData")
data <- as.data.frame(data)
load("~/Dropbox (Personal)/Bioinformatics tools/ICR genes/ICR_genes.RData")

# Sidra-LUMC
load("./Processed_Data/RNASeq/Trimmed_p/Raw/JSREP.Complete.gene.filtered.Rdata")

# Analysis

# drop rows(genes) without approved gene-symbol in TCGA dataset
RNAseq.genes = rownames(data)
info.genes = rownames(geneInfo)
available.genes = unique(RNAseq.genes[which(RNAseq.genes %in% info.genes)])
geneInfo = geneInfo[which(rownames(geneInfo) %in% available.genes),]                                                      # drop the genes without RNAseq.DATA
data = data[which(rownames(data) %in% available.genes),]         # drop the genes without info

# Make gene order identical in both datasets
intersected_genes = intersect(rownames(data), rownames(expression.filtered))
ICR_genes %in% intersected_genes

data = data[intersected_genes,] # TCGA-COAD
expression.filtered = expression.filtered[intersected_genes,] # Sidra-LUMC
geneInfo = geneInfo[intersected_genes,]  # GeneInfo

dim(data)
dim(expression.filtered)
dim(geneInfo)

# cbind both matrices
RNAseqData_Merged = cbind(data, expression.filtered)
RNAseqData_Merged = as.matrix(RNAseqData_Merged)
mode(RNAseqData_Merged) = "numeric"
sum(is.na(RNAseqData_Merged))   # check if there are no NA's : 0

# Normalization
RNASeq.expr.set = newSeqExpressionSet(RNAseqData_Merged, featureData = geneInfo)                              # Create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "gcContent"] = as.numeric(geneInfo[, "gcContent"])                                   # Make sure gcContenet is numeric
RNASeq.expr.set = withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE)       # Removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set = betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)                   # Removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM = log(RNAseqData_Merged + .1) + offst(RNASeq.expr.set)                                            # Apply the Edaseq Offset
RNASeq.NORM = floor(exp(RNASeq.NORM) - .1)                                                                    # Return non decimal values

#Quantile normalisation RNA
RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)                                         # Quantile normalize
RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                                             # Return non decimal values
rownames(RNASeq.NORM.quantiles) <- rownames(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)

dim(RNASeq.NORM.quantiles)  # Check number of rows and columns: 18291, 911

# Save complete file
dir.create("./Processed_Data/RNASeq/030_Merged_TCGA_COAD_Sidra_LUMC", showWarnings = FALSE)
save(RNASeq.NORM.quantiles, geneInfo, file = "./Processed_Data/RNASeq/030_Merged_TCGA_COAD_Sidra_LUMC/TCGA_COAD_Sidra_LUMC_All_Samples_EDASEQ_QN.Rdata")

# Filter samples
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
primary_Sidra_LUMC = colnames(RNASeq.QN.LOG2)
load("../NGS_Data_TCGA_COAD_Jessica/Analysis/Deconvolution_and_GSEA/Bindea_ORIG_ES.Rdata")
primary_TCGA_COAD = colnames(ES)

All_primary_samples = c(primary_Sidra_LUMC, primary_TCGA_COAD)  # length(All_primary_samples) = 787

RNASeq.NORM.quantiles = RNASeq.NORM.quantiles[, All_primary_samples]
dim(RNASeq.NORM.quantiles) # 18291   787

# log transformation
RNASeq.QN.LOG2=log(RNASeq.NORM.quantiles+1,2)

save(RNASeq.NORM.quantiles, RNASeq.QN.LOG2, geneInfo, file = "./Processed_Data/RNASeq/030_Merged_TCGA_COAD_Sidra_LUMC/TCGA_COAD_Sidra_LUMC_Primary_tumor_EDASEQ_QN_LOG2.Rdata")
