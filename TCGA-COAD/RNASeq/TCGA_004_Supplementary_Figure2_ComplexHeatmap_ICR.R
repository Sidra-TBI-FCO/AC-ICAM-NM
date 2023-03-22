#################################################################
###
### This script creates heatmaps for ICR genes using EDASeq 
### normalized RNASeq data. Samples are ordered by ICR 
### score and a heatmap is created for ICR genes including 
### ICR cluster allocation for k = 3, 4 and 6 and the proposed
### HML classfication.
### Heatmaps are saved as png files at location: 
### 
###
#################################################################

## Create Heatmap of TCGA RNASeq ICR genes
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
download.method = "Biolinks"

# Load data
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))

# Create folders
dir.create("./Figures/",showWarnings = FALSE)
dir.create("./Figures/Heatmaps", showWarnings = FALSE)
dir.create("./Figures/Heatmaps/ICR Heatmaps", showWarnings = FALSE)

if(download.method == "TCGA_Assembler"){
  load("./Analysis/ICR Consensus Clustering/COAD_ICR_cluster_assignment_k2-6.Rdata")
  load("./Processed_Data/RNASeq/TCGA_Assembler_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
  load("./Analysis/016_CMS_Classification/TCGA_Assembler_Rfcms.Rdata")
}

if(download.method == "Biolinks"){
  load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
  load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
  load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
}

load("~/Dropbox (SMCP)/TCGA pancancer-germline/New Master Data/Master.file.new.RS.v3.Rdata")

# Analysis
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)

ICR_subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% ICR_genes, ])
ICR_subset_RNAseq_log2 = ICR_subset_RNAseq_log2[,ICR_genes]

table_cluster_assignment$ICR_HML = factor(table_cluster_assignment$ICR_HML, levels = c("ICR High", "ICR Medium", "ICR Low"))

#table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICRscore),]
table_cluster_assignment = table_cluster_assignment[order(table_cluster_assignment$ICR_HML, decreasing = TRUE),]
table_cluster_assignment$Sample_ID = rownames(table_cluster_assignment)
table_cluster_assignment$CMS = Rfcms$RF.predictedCMS[match(table_cluster_assignment$Sample_ID, rownames(Rfcms))]
table_cluster_assignment$CMS[which(is.na(table_cluster_assignment$CMS))] = "mixed"
colnames(table_cluster_assignment)[which(colnames(table_cluster_assignment) == "HML_cluster")] = "ICR_HML"

# MANTIS
pheno.data.COAD = pheno.data[which(pheno.data$Patient.ID %in% substring(table_cluster_assignment$Sample_ID, 1, 12)),]
MANTIS = pheno.data.COAD[, c("Patient.ID", "MANTIS_score")]
colnames(MANTIS) = c("Patient_ID", "MANTIS.score")
MANTIS$MSI[which(MANTIS$MANTIS.score > 0.4)] = "MSI-H"
MANTIS$MSI[which(MANTIS$MANTIS.score <= 0.4)] = "MSS"

dir.create("./Processed_Data/External_Data", showWarnings = FALSE)
dir.create("./Processed_Data/External_Data/From_TCGA_Masterfile", showWarnings = FALSE)

save(MANTIS, file = paste0("./Processed_Data/External_Data/From_TCGA_Masterfile/", download.method, "_MANTIS.Rdata"))

table_cluster_assignment$MSI = MANTIS$MSI[match(substring(rownames(table_cluster_assignment), 1, 12),
                                                MANTIS$Patient_ID)]
table_cluster_assignment$MSI[which(is.na(table_cluster_assignment$MSI))] = "Not determined"

table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(RNASeq.QN.LOG2)),]
sample_order = rownames(table_cluster_assignment)

Expression.matrix = t(ICR_subset_RNAseq_log2[sample_order,])

# z-score Expression.matrix
Expression.matrix.z = Expression.matrix
for(j in 1: nrow(Expression.matrix.z))  {
  Expression.matrix.z[j,] = (Expression.matrix[j,]-mean(Expression.matrix[j,]))/sd(Expression.matrix[j,]) # z-score the enrichment matrix
}

table_cluster_assignment$ICR_HML = factor(table_cluster_assignment$ICR_HML, levels = c("ICR High", "ICR Medium", "ICR Low"))
table_cluster_assignment$CMS = factor(table_cluster_assignment$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
table_cluster_assignment$MSI = factor(table_cluster_assignment$MSI, levels = c("MSI-H", "MSS", "Not determined"))


ha = HeatmapAnnotation(`ICR cluster` = table_cluster_assignment$ICR_HML,
                       `CMS` = table_cluster_assignment$CMS,
                       MSI = table_cluster_assignment$MSI,
                       col = list(`ICR cluster` = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "white", "Not determined" = "lightgrey"),
                                  `CMS` = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                                    "mixed" = "lightgrey")))

png(filename = paste0("./Figures/Heatmaps/ICR Heatmaps/", download.method ,"_ICR_ComplexHeatmap.png"), res = 600,
    width = 5.5, height = 3.5, units = "in")
Heatmap(Expression.matrix.z, cluster_rows = FALSE, cluster_columns = FALSE,
        show_column_names = FALSE, top_annotation = ha, name = "Expression\n z score"
)
dev.off()
