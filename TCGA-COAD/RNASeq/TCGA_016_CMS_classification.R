#################################################################
###
### CMS classification of colon samples
###
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

library(devtools)
install_github("Sage-Bionetworks/CMSclassifier")
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
install.packages("synapser")

#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")

library(synapser)
library(CMSclassifier)

# Set parameters
download.method = "Biolinks"

# Load data
if(download.method == "TCGA_Assembler"){
  load("./Processed_Data/RNASeq/TCGA_Assembler_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata")
}
if(download.method == "Biolinks"){
  load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
}

load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/GeneInfo/geneInfo.July2017.RData")
CMS_genes = CMSclassifier::listModelGenes("RF")
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)

#Sample----------273 genes used for CMS classification, these need to be included to make a CMS prediction
#load("./Analysis/CMS/CMS_SampleData.Rdata") #comment out later
#sampleDataSubset = sampleData[,colnames(sampleData) %in% CMS_genes]
#Test_Rfcms = CMSclassifier::classifyCMS(t(sampleDataSubset),method="RF")[[3]]
#---this works
#sampleDataSubset_withNAs = sampleDataSubset
#sampleDataSubset_withNAs$`4628` = NA
#sampleDataSubset_withNAs$`123920` = NA
#Test2_Rfcms = CMSclassifier::classifyCMS(t(sampleDataSubset_withNAs),method="RF")[[3]]
#----------------this works too

RNASeq.data.entrez = RNASeq.QN.LOG2
row.names(RNASeq.data.entrez) = geneInfo$EntrezID[match(row.names(RNASeq.data.entrez), row.names(geneInfo))]
RNASeq.data.entrez = as.data.frame(t(RNASeq.data.entrez))

available_CMS_genes = CMS_genes[which(CMS_genes %in% colnames(RNASeq.data.entrez))]
unavailable_CMS_genes = CMS_genes[-which(CMS_genes %in% colnames(RNASeq.data.entrez))]

# Adding unavailable CMS genes to expression data (as NA's)
RNASeq.data.entrez$`8857` = NA
RNASeq.data.entrez$`53826` = NA
RNASeq.data.entrez$`6935` = NA
RNASeq.data.entrez$`23475` = NA
RNASeq.data.entrez$`4360` = NA

#RNASeq.data.entrez$`53826` = NA
#RNASeq.data.entrez$`196051` = NA
#RNASeq.data.entrez$`57608` = NA
#RNASeq.data.entrez$`6935` = NA
#RNASeq.data.entrez$`1687` = NA
#RNASeq.data.entrez$`23475` = NA
#RNASeq.data.entrez$`4360` = NA
#RNASeq.data.entrez$`23111` = NA
#RNASeq.data.entrez$`284119` = NA
#RNASeq.data.entrez$`2983` = NA
#RNASeq.data.entrez$`2982` = NA

RNASeq.data.entrez = RNASeq.data.entrez[,-which(is.na(colnames(RNASeq.data.entrez)))]
RNASeq.data.entrez = RNASeq.data.entrez[,which(colnames(RNASeq.data.entrez) %in% CMS_genes)]

Rfcms = CMSclassifier::classifyCMS(t(RNASeq.data.entrez),method="RF")[[3]]
SScms = CMSclassifier::classifyCMS(t(RNASeq.data.entrez),method="SSP")[[3]]
table(Rfcms$RF.predictedCMS, exclude = NULL)
prop.table(table(Rfcms$RF.predictedCMS, exclude = NULL))*100
table(SScms$SSP.predictedCMS, exclude = NULL)
prop.table(table(SScms$SSP.predictedCMS, exclude = NULL))*100

dir.create("./Analysis/016_CMS_Classification", showWarnings = FALSE)
save(Rfcms, SScms, file = paste0("./Analysis/016_CMS_Classification/", download.method, "_Rfcms.Rdata"))

