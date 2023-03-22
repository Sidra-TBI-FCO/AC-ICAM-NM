####################################################################
###
### This Script clusters the RNASeq data that have been
### normalized using EDASeq. Additionally, optimal kalinsky
### is calculated.
### 
### Input data:
### ("~/Dropbox (TBI-Lab)/DB-LAB/Projects/JSREP Colon Cancer/Data/Complete Cohort/At Sidra/RNASeq Data/Normalized/JSREP.dataset.EDAseq.QN.HPC.Rdata")
### Output data are saved as Rdata file:
### ("~/Dropbox (TBI-Lab)/DB-LAB/Projects/JSREP Colon Cancer/Data/Complete Cohort/At Sidra/RNASeq Data/ICR Consensus Clustering/JSREP.dataset.ICR.cluster.assignment.k2-6.Rdata"
### which includes: (1) table_cluster_assignment and (2) optimal.calinsky.
#####################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/stefanofunctions.R"))                                                                 # Used for calinsky function and plot

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper")
required.bioconductor.packages = c("ConsensusClusterPlus", "clue")
ipak(required.packages)
ibiopak(required.bioconductor.packages)

# Set Parameters
set.seed(6)

# Load data
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/ICR genes/ICR_genes.RData")
clinical_data = read.csv("./Processed_Data/External_Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv", stringsAsFactors = FALSE)

# Create folders
dir.create("./Analysis",showWarnings = FALSE)
dir.create("./Analysis/ICR Consensus Clustering",showWarnings = FALSE)

# Log transform
RNASeq.QN.LOG2 = log(filtered.norm.RNAseqData +1, 2)

## Perform clusting
ICR_subset_RNAseq_log2 = t(RNASeq.QN.LOG2[row.names(RNASeq.QN.LOG2) %in% ICR_genes, ])
ddist = dist(ICR_subset_RNAseq_log2)

setwd("./Analysis/ICR Consensus Clustering")

ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 6,                                                                              # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                             # set repeats
                                               title=paste0("TCGA.GDC_Biolinks.ICR.reps5000"),                                       # Output filename (no path)
                                               clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE)

outputfiles = list.files("TCGA.GDC_Biolinks.ICR.reps5000", full.names = TRUE)
class_files = outputfiles[grep("consensusClass", outputfiles)]

N.files = length(class_files)
table_cluster_assignment = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))

for (i in 1:N.files){
  file = paste0("./", class_files[i])
  consensus_class = read.csv(file = file,header=FALSE)
  group = paste0("Group_k",i+1)
  colnames(consensus_class) = c("PatientID", group)
  rownames(consensus_class) = consensus_class$PatientID
  consensus_class$PatientID = NULL
  table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))]
  
  transl_table_ICR_cluster = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)
  colnames(transl_table_ICR_cluster) = c(group,"mean_ICRscore")
  transl_table_ICR_cluster = cbind(transl_table_ICR_cluster[order(transl_table_ICR_cluster$mean_ICRscore),],ICR_name=paste0("ICR",c(1:(i+1))))
  
  ICR_cluster = paste0("ICR_cluster_k",i+1)
  table_cluster_assignment[,ICR_cluster] = transl_table_ICR_cluster$ICR_name[match(table_cluster_assignment[,group],
                                                                                   transl_table_ICR_cluster[,group])]
}

#calinsky
sHc <- hclust(ddist, method = "ward.D2")
aCalinsky <- calinsky(sHc,gMax=10)
pdf(file = "./TCGA_GDC_Biolinks_ICR_cluster_assignment_k2-6.Calinsky.pdf", width = 16, height = 6)
plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
dev.off()
optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))

#save data
outputname = paste0("./TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")

save(table_cluster_assignment,optimal.calinsky, file = outputname)


# Chose optimal k = 3 and add ICR Low, ICR Medium and ICR High annotation
table_cluster_assignment$ICR_HML = NA
table_cluster_assignment$ICR_HML[which(table_cluster_assignment$ICR_cluster_k3 == "ICR1")] = "ICR Low"
table_cluster_assignment$ICR_HML[which(table_cluster_assignment$ICR_cluster_k3 == "ICR2")] = "ICR Medium"
table_cluster_assignment$ICR_HML[which(table_cluster_assignment$ICR_cluster_k3 == "ICR3")] = "ICR High"
table_cluster_assignment$ICR_HML = factor(table_cluster_assignment$ICR_HML, levels = c("ICR Low", "ICR Medium", "ICR High"))

save(table_cluster_assignment, optimal.calinsky, file = outputname)
