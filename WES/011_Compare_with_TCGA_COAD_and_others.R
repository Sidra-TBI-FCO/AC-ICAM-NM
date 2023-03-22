
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr")
ipak(required.packages)

# Load data
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
Sidra_LUMC = read.csv("./Analysis/WES/001_Oncoplot/Gene_Count_MAF.csv", stringsAsFactors = FALSE)
TCGA_COAD = read.csv("../NGS_Data_TCGA_COAD_Jessica/Analysis/WES/001_Oncoplot/Gene_Count_MAF.csv",
                     stringsAsFactors = FALSE)
Giannakis_nonhypermutated = read.csv("./Processed_Data/External/Gene_collections/Giannakis_MutSig_analysis_nonhypermutated.csv", stringsAsFactors = FALSE)
Giannakis_hypermutated = read.csv("./Processed_Data/External/Gene_collections/Giannakis_MutSig_analysis_hypermutated.csv", stringsAsFactors = FALSE)
load("./Analysis/WES/011_Compare_with_TCGA/NHS_HPFS_Data_Frequency_Nonsilent_mutations_in_Colon_only_482_patients.Rdata")
Grasso = read.csv("./Processed_Data/External/Gene_collections/Grasso_Supplementary_Table6.csv", stringsAsFactors = FALSE)
COSMIC = read.csv("./Processed_Data/External/Gene_collections/COSMIC_Colon_Cancer_Adenocarcinoma_Genes_with_mutation_GRCh37_22_Sept_2020.csv", stringsAsFactors = FALSE)

QGP_file = read.csv("./Processed_Data/External/Gene_collections/QGPC_Master_Gene_List_1278_RELEASE1_test1._MAY7_ED.csv", stringsAsFactors = FALSE)

Bailey_genes = read.csv("./Processed_Data/External/Gene_collections/Bailey_Cell_S2.csv", stringsAsFactors = FALSE) 
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
Oncogenic_mediators_Colaprico = read.csv("./Processed_Data/External/Gene_collections/Supplementary_table_6_Colaprico.csv", stringsAsFactors = FALSE)

# Prepare COSMIC data
COSMIC = COSMIC[-grep("_ENST", COSMIC$Gene.name),]
COSMIC$Frequency = COSMIC$Mutated.samples / COSMIC$Samples.tested

# Split for hypermutated and non-hypermutated
frequency_df$Patient_ID = as.character(frequency_df$Patient_ID)
hypermutated = frequency_df$Patient_ID[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)]
nonhypermutated = frequency_df$Patient_ID[-which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)]

finalMafFiltered$Tumor_Sample_Barcode = as.character(finalMafFiltered$Tumor_Sample_Barcode)
h_MAF = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% hypermutated),]
nh_MAF = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% nonhypermutated),]

# Gene count ranking
# Hypermutated
h_MAF_gene_Count= h_MAF %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count)) 

h_MAF_gene_Count$Percentage = round((h_MAF_gene_Count$count / length(hypermutated)) * 100, 1)
write.csv(h_MAF_gene_Count, file = "./Analysis/WES/001_Oncoplot/Gene_Count_MAF_hypermutated.csv")

# Nonhypermutated
nh_MAF_gene_Count= nh_MAF %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count)) 

nh_MAF_gene_Count$Percentage = round((nh_MAF_gene_Count$count / length(nonhypermutated)) * 100, 1)
write.csv(nh_MAF_gene_Count, file = "./Analysis/WES/001_Oncoplot/Gene_Count_MAF_nonhypermutated.csv")

# Add to overview file
Sidra_LUMC = Sidra_LUMC[, c("Hugo_Symbol", "Percentage")]
colnames(Sidra_LUMC) = c("Hugo_Symbol", "Percentage_in_Sidra_LUMC_cohort")
Sidra_LUMC$Percentage_in_nonhypermutated_Sidra_LUMC_cohort = nh_MAF_gene_Count$Percentage[match(Sidra_LUMC$Hugo_Symbol, nh_MAF_gene_Count$Hugo_Symbol)]
Sidra_LUMC$Percentage_in_hypermutated_Sidra_LUMC_cohort = h_MAF_gene_Count$Percentage[match(Sidra_LUMC$Hugo_Symbol, h_MAF_gene_Count$Hugo_Symbol)]

# Compare with TCGA COAD
Sidra_LUMC$Percentage_in_TCGA_COAD = NA
Sidra_LUMC$Percentage_in_TCGA_COAD = TCGA_COAD$Percentage[match(Sidra_LUMC$Hugo_Symbol, TCGA_COAD$Hugo_Symbol)]

# COSMIC
Sidra_LUMC$Percentage_in_COSMIC_Colon_Adenocarcinoma = round(COSMIC$Frequency[match(Sidra_LUMC$Hugo_Symbol,
                                                                                    COSMIC$Gene.name)] *100, 1)

# With Giannakis (Hypermutated tumors are defined as those with a mutation rate > 12/Mb)
Giannakis_nonhypermutated$Percentage = round((Giannakis_nonhypermutated$npat / 488) *100, 1)
Sidra_LUMC$Giannakis_MutSig_Percentage_in_NHS_HPFS_nonhypermutated = Giannakis_nonhypermutated$Percentage[match(Sidra_LUMC$Hugo_Symbol,
                                                                                                                Giannakis_nonhypermutated$Gene)]
Giannakis_hypermutated$Percentage = round((Giannakis_hypermutated$npat / 131) *100, 1)
Sidra_LUMC$Giannakis_MutSig_Percentage_in_NHS_HPFS_hypermutated = Giannakis_hypermutated$Percentage[match(Sidra_LUMC$Hugo_Symbol,
                                                                                                          Giannakis_hypermutated$Gene)]

# With Giannakis (colon only)
df$Frequency = round(df$Frequency *100, 1)
Sidra_LUMC$NHS_PFHS = df$Frequency[match(Sidra_LUMC$Hugo_Symbol, df$Gene)]

# With Grasso (MSI-H subgroup and MSS subgroup)
Grasso_MSI_H_genes = Grasso$Gene[which(Grasso$TCGA.NHS.HPFS.CRC.MSI.high.MutSigCV.q.0.01 == "Yes")]
Sidra_LUMC$Grasso_MutSig_in_MSI_H_group = NA
Sidra_LUMC$Grasso_MutSig_in_MSI_H_group[which(Sidra_LUMC$Hugo_Symbol %in% Grasso_MSI_H_genes)] = "yes"

Grasso_MSS_genes = Grasso$Gene[which(Grasso$TCGA.NHS.HPFS.CRC.MSS.MutSigCV.q.0.01 == "Yes")]
Sidra_LUMC$Grasso_MutSig_in_MSS_group = NA
Sidra_LUMC$Grasso_MutSig_in_MSS_group[which(Sidra_LUMC$Hugo_Symbol %in% Grasso_MSS_genes)] = "yes"

#MSK_IMPACT_genes = QGP_file$Source.Gene.Symbol[which(QGP_file$MSK.IMPACT_OncoKB == "Yes")]
#Sidra_LUMC$MSK_IMPACT_test = NA
#Sidra_LUMC$MSK_IMPACT_test[which(Sidra_LUMC$Hugo_Symbol %in% MSK_IMPACT_genes)] = "yes"

# Identified drivers
PANCAN_genes = Bailey_genes$Gene[which(Bailey_genes$Cancer.type == "PANCAN")]
Sidra_LUMC$Gene_in_Bailey_PANCAN = "no"
Sidra_LUMC$Gene_in_Bailey_PANCAN[which(Sidra_LUMC$Hugo_Symbol %in% PANCAN_genes)] = "yes"

COADREAD_genes = Bailey_genes$Gene[which(Bailey_genes$Cancer.type == "COADREAD")]
Sidra_LUMC$Gene_in_Bailey_COADREAD = "no"
Sidra_LUMC$Gene_in_Bailey_COADREAD[which(Sidra_LUMC$Hugo_Symbol %in% COADREAD_genes)] = "yes"

Sidra_LUMC$Gene_in_QGP_collection = "no"
Sidra_LUMC$Gene_in_QGP_collection[which(Sidra_LUMC$Hugo_Symbol %in% QGPC_genes)] = "yes"

QGP_file = QGP_file[which(QGP_file$Source.Gene.Symbol %in% QGPC_genes),]
QGP_file$Source = NA
rownames(QGP_file) = QGP_file$Source.Gene.Symbol
QGP_file$QGP_CANCER_RELATED_LIST..28.ACMG.NON.CANCER.and.24.CHIARA.NON.CANCER.EXCLUDED._1226 = NULL

i = 1
for (i in 1:nrow(QGP_file)){
  gene = rownames(QGP_file)[i]
  sources = colnames(QGP_file)[which(QGP_file[gene,] %in% c("YES", "Yes"))]
  QGP_file[gene, "Source"] = paste(sources, collapse = " & ")
}

Sidra_LUMC$Source_QGP_collection = QGP_file$Source[match(Sidra_LUMC$Hugo_Symbol, QGP_file$Source.Gene.Symbol)]


COAD_genes = Oncogenic_mediators_Colaprico$Gene[which(Oncogenic_mediators_Colaprico$Cancertype == "COAD")]
Sidra_LUMC$Gene_in_Colaprico_COAD = "no"
Sidra_LUMC$Gene_in_Colaprico_COAD[which(Sidra_LUMC$Hugo_Symbol %in% COAD_genes)] = "yes"

Sidra_LUMC$Artifact_gene = "no"
artifacts = c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG",paste0("OR",1:9))
Sidra_LUMC$Artifact_gene[which(substring(Sidra_LUMC$Hugo_Symbol, 1, 3) %in% artifacts)] = "yes"
artifacts2 = c("PLIN","CELA","SRA1")
Sidra_LUMC$Artifact_gene[which(substring(Sidra_LUMC$Hugo_Symbol, 1, 4) %in% artifacts2)] = "yes"
artifacts_genes = c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")
Sidra_LUMC$Artifact_gene[which(Sidra_LUMC$Hugo_Symbol %in% artifacts_genes)] = "yes"

dir.create("./Analysis/WES/011_Compare_with_TCGA", showWarnings = FALSE)
write.csv(Sidra_LUMC, file = "./Analysis/WES/011_Compare_with_TCGA/Oct_2021_Sidra_LUMC_frequency_compared_with_TCGA_NHS_PFHS.csv",
          row.names = TRUE)
