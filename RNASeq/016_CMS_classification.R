#################################################################
###
### CMS classification of colon samples
###
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

library(devtools)
install_github("Sage-Bionetworks/CMSclassifier")
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
install.packages("synapser")

#source('http://depot.sagebase.org/CRAN.R')
#pkgInstall("synapseClient")

library(synapser)
library(CMSclassifier)

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/GeneInfo/geneInfo.July2017.RData")
CMS_genes = CMSclassifier::listModelGenes("RF")

RNASeq.data.entrez = RNASeq.QN.LOG2
row.names(RNASeq.data.entrez) = geneInfo$EntrezID[match(row.names(RNASeq.data.entrez), row.names(geneInfo))]
RNASeq.data.entrez = as.data.frame(t(RNASeq.data.entrez))

available_CMS_genes = CMS_genes[which(CMS_genes %in% colnames(RNASeq.data.entrez))]
unavailable_CMS_genes = CMS_genes[-which(CMS_genes %in% colnames(RNASeq.data.entrez))]

# Adding unavailable CMS genes to expression data (as NA's)
RNASeq.data.entrez$`53826` = NA
RNASeq.data.entrez$`196051` = NA
RNASeq.data.entrez$`57608` = NA
RNASeq.data.entrez$`6935` = NA
RNASeq.data.entrez$`1687` = NA
RNASeq.data.entrez$`23475` = NA
RNASeq.data.entrez$`4360` = NA
RNASeq.data.entrez$`23111` = NA
RNASeq.data.entrez$`284119` = NA
RNASeq.data.entrez$`2983` = NA
RNASeq.data.entrez$`2982` = NA


RNASeq.data.entrez = RNASeq.data.entrez[,-which(is.na(colnames(RNASeq.data.entrez)))]
RNASeq.data.entrez = RNASeq.data.entrez[,which(colnames(RNASeq.data.entrez) %in% CMS_genes)]

Rfcms = CMSclassifier::classifyCMS(t(RNASeq.data.entrez),method="RF")[[3]]
SScms = CMSclassifier::classifyCMS(t(RNASeq.data.entrez),method="SSP")[[3]]

table(Rfcms$RF.predictedCMS, exclude = NULL)
prop.table(table(Rfcms$RF.predictedCMS, exclude = NULL))*100
table(SScms$SSP.predictedCMS, exclude = NULL)
prop.table(table(SScms$SSP.predictedCMS, exclude = NULL))*100

dir.create("./Analysis/Trimmed_p/016_CMS_Classification", showWarnings = FALSE)
save(Rfcms, SScms, file = "./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

