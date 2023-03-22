
#script for processing and normalistaton of RNAseq data from subreads
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

#read text file
input.file = "./RNASeq_HPC/Create expression matrix/Trimmed_p/counts.txt"
subread.command = readLines (input.file,n = 1)
raw.data = read.csv(input.file,stringsAsFactors = FALSE,skip = 1,sep = "\t")
  
#extract relevant data
gene.data = raw.data[,c(1:6)]
expression.data = raw.data[,-c(2:6)]

#expression.data = expression.data[,-ncol(expression.data)]
rownames(expression.data) = expression.data$Geneid
expression.data$Geneid = NULL
colnames(expression.data) = gsub(x = colnames(expression.data),pattern = "...TRIMMED.BAM.RNASeq_COAD_LUMC_SIDRA_",replacement = "")
colnames(expression.data) = gsub(x = colnames(expression.data),pattern = ".bam",replacement = "")
colnames(expression.data) = gsub(x = colnames(expression.data), pattern = "1LM1", replacement = "LM1")
colnames(expression.data) = gsub(x = colnames(expression.data), pattern = "2LM2", replacement = "LM2")
colnames(expression.data) = gsub(x = colnames(expression.data), pattern = ".trimmed", replacement = "")
colnames(expression.data) = substring(colnames(expression.data), 1, 6)
colnames(expression.data) = gsub(x = colnames(expression.data), pattern = "\\.", replacement = "_")
expression.matrix = as.matrix(expression.data)
#save as R datafile
save (gene.data,subread.command,expression.matrix,file="./Data/RNASeq/Trimmed_p/Raw/JSREP.Complete.dataset.HPC.Rdata")

# normalizationa
# dependencies
ibiopak("EDASeq")
ibiopak("base64enc")
ibiopak("preprocessCore")
load ("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/GeneInfo/geneInfo.Sept2018.RData")                                                               

# cleanup expression data
geneInfo = as.data.frame(geneInfo)
geneInfo = geneInfo[!is.na(geneInfo[,1]),]                                                                          # drop the genes without information and convert to data frame
available.genes = unique(rownames(expression.matrix)[which(rownames(expression.matrix) %in% rownames(geneInfo))])   # genes in geneinfo and expression.matrix
geneInfo = geneInfo[available.genes,]                                                                               # drop the GCcontent info on unavailable genes
expression.filtered = expression.matrix[available.genes,]                                                           # drop the genes without information from the RNAseq.DATA
mode(expression.filtered) = "numeric"
dim(expression.filtered)  #rows = 19460, columns =390

save(expression.filtered, file = "./Data/RNASeq/Trimmed_p/Raw/JSREP.Complete.gene.filtered.Rdata")

# EDAseq 
RNASeq.expr.set = newSeqExpressionSet(expression.filtered , featureData = geneInfo)                     # create a new SeqExpressionSet object.
fData(RNASeq.expr.set)[, "gcContent"] = as.numeric(geneInfo[, "gcContent"])                             # make sure gcContenet is numeric
RNASeq.expr.set = withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE) # removes lane gene specific effects, for example effects related to gene length or GC content
RNASeq.expr.set = betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)             # removes effect related to in between lane distributional differences, as sequencing depth
RNASeq.NORM =  log(expression.filtered + .1) + offst(RNASeq.expr.set)                                   # apply the Edaseq Ofset
RNASeq.NORM =  floor(exp(RNASeq.NORM) - .1)                                                             # return non decimal values

# Quantile normalisation RNA
RNASeq.NORM.quantiles = normalize.quantiles(RNASeq.NORM)                                                # Quantile normalize
RNASeq.NORM.quantiles = floor(RNASeq.NORM.quantiles)                                                    # return non decimal values
rownames(RNASeq.NORM.quantiles) = rownames(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) = colnames(RNASeq.NORM)

# log transformation
RNASeq.QN.LOG2=log(RNASeq.NORM.quantiles+1,2)
dir.create("./Data/RNASeq/Trimmed_p/Normalized", showWarnings = FALSE)
save (RNASeq.QN.LOG2,RNASeq.NORM.quantiles,geneInfo,file="./Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")
