
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ComplexHeatmap")
ipak(required.packages)

# Load data
load("./Analysis/WES/008_matrix_for_Oncoplot/matrix_selected_for_oncoplot.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/External/Gene_collections/QGPC_genes.Rdata")
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_thresholded_by_genes.Rdata")

# Set parameters
artifacts_removed = "artifacts_removed"
QGP_only = "QGP_only"

# Analysis
mat[1:3, 1:3]

if(artifacts_removed == "artifacts_removed"){
  artifacts = c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG",paste0("OR",1:9))
  mat = mat[-which(substring(rownames(mat), 1, 3) %in% artifacts),]
  artifacts2 = c("PLIN","CELA","SRA1")
  if(substring(rownames(mat), 1, 4) %in% artifacts >0){
    mat = mat[-which(substring(rownames(mat), 1, 4) %in% artifacts),]
    }
  artifacts_genes = c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")
  mat = mat[-which(rownames(mat) %in% artifacts_genes),]
}
if(QGP_only == "QGP_only"){
  mat = mat[which(rownames(mat) %in% QGPC_genes),]
}

mat = mat[c(1:21, which(rownames(mat) == "POLE")),]  # For paper: 1:21

df_GISTIC = df_GISTIC[rownames(mat),colnames(mat)]

i=1
j=1
for(i in 1:nrow(mat)){
  row = rownames(mat)[i]
  for (j in 1:ncol(mat)){
    col = colnames(mat)[j]
    CNV = df_GISTIC[row, col]
    CNV = as.numeric(CNV)
    if(CNV == 2 & mat[row,col] == "MUT"){
      print(paste0("MUT", row, col))
      mat[row,col] = paste0(mat[row, col], ";AMP")
    }
    if(CNV == -2){
      mat[row, col] = paste0(mat[row, col], "HOMDEL")
    }
    if(CNV == 2 & mat[row,col] == ""){
      mat[row,col] = paste0(mat[row, col], "AMP")
    }
    
  }
}

# Define functions for OncoPrint
col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "black")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#e6e6e6", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)

column_title = ""  #"Sidra-LUMC cohort \n Most frequently mutated genes"
heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "MUT"), 
                            labels = c("Deep deletion", "Amplification", "Mutation"))

# Prepare for annotation
MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
Rfcms = Rfcms[which(substring(rownames(Rfcms), 1, 4) %in% colnames(mat)),]
Rfcms$RF.predictedCMS = as.character(Rfcms$RF.predictedCMS)
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"

df = data.frame(Sample_ID = substring(rownames(Rfcms), 1, 4), CMS = Rfcms$RF.predictedCMS, 
                ICR = NA, MSI = NA, Mutational_load = NA)
rownames(df) = df$Sample_ID
df$MSI = MANTIS$MSI[match(df$Sample_ID, MANTIS$Sample_ID)]
df$Mutational_load = frequency_df$Non_silent_Mutation_frequency[match(substring(df$Sample_ID, 1, 3), frequency_df$Patient_ID)]
df$ICR = table_cluster_assignment$ICR_HML[match(df$Sample_ID, substring(rownames(table_cluster_assignment), 1, 4))]

df$ICR = factor(df$ICR, levels = c("ICR High", "ICR Medium", "ICR Low"))
df$CMS = factor(df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
df$MSI = factor(df$MSI, levels = c("MSI-H", "MSS"))
  
# Define order based on mutational load
df = df[order(df$Mutational_load),]
sample_order = rownames(df)
mat = mat[,sample_order]

# Heatmap annotation
ha = HeatmapAnnotation(`Mutational load` = anno_barplot(df$Mutational_load, baseline = 0),
                       ICR = df$ICR,
                       CMS = df$CMS,
                       MSI = df$MSI,
                       col = list(ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  CMS = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                          "mixed" = "lightgrey"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "white"))
)

# Plotting
dir.create("./Figures/WES/009_Oncoprint", showWarnings = FALSE)
#png(paste0("./Figures/WES/009_Oncoprint/Aug_2020_OncoPrint_POLE_", artifacts_removed,".png"), res = 600, width = 11, height = 7.5,
 #   units = "in")
pdf(paste0("./Figures/WES/009_Oncoprint/Fig3c_Nov_2020_OncoPrint_CNV_POLE_", artifacts_removed, QGP_only, ".pdf"), width = 12, height = 5)
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          top_annotation = ha,
          row_order = rownames(mat),
          column_order = sample_order)
dev.off()

genes_Sidra_LUMC_Oncoprint = rownames(mat)
dir.create("./Analysis/WES/009_Oncoprint", showWarnings = FALSE)
save(genes_Sidra_LUMC_Oncoprint, 
     file = paste0("./Analysis/WES/009_Oncoprint/Genes_Oncoprint_", artifacts_removed, QGP_only, ".Rdata"))
