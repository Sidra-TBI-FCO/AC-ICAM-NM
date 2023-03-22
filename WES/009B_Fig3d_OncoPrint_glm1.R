
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ggplot2", "ComplexHeatmap")
ipak(required.packages)

# Set parameters
version = "glm_1" # "glm_1"
subset = "all_patients" # "hypermutated" "all_patients"

# Load data
load(paste0("./Analysis/WES/008_matrix_for_Oncoplot/008", version, "_matrix_for_oncoplot.Rdata"))
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/WES/GISTIC/001_Data_prep/GISTIC_281_patients_thresholded_by_genes.Rdata")

# Analysis
#df_GISTIC = df_GISTIC[rownames(mat),colnames(mat)]
# i=1
#j=1
#for(i in 1:nrow(mat)){
 # row = rownames(mat)[i]
  #for (j in 1:ncol(mat)){
   # col = colnames(mat)[j]
   # CNV = df_GISTIC[row, col]
   # CNV = as.numeric(CNV)
   # if(CNV == 2 & mat[row,col] == "MUT"){
   #   print(paste0("MUT", row, col))
   #   mat[row,col] = paste0(mat[row, col], ";AMP")
   # }
   # if(CNV == -2 & mat[row,col] == "MUT"){
   #   mat[row, col] = paste0(mat[row, col], ";HOMDEL")
  #  }
  #  if(CNV == -2 & mat[row,col] == ""){
   #   mat[row, col] = paste0(mat[row, col], "HOMDEL")
   # }
   # if(CNV == 2 & mat[row,col] == ""){
   #   mat[row,col] = paste0(mat[row, col], "AMP")
   # }
    
  #}
#}



# Define functions for OncoPrint
col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "black")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#e6e6e6", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)

#column_title = "Sidra-LUMC cohort \n BRCA1, BRCA2, FA, NON-BRCA HR somatic mutations"
column_title = ""
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
df$Mutational_load_per_MB = frequency_df$Nonsilent_mutational_burden_per_Mb[match(substring(df$Sample_ID, 1, 3), frequency_df$Patient_ID)]
df$ICR = table_cluster_assignment$ICR_HML[match(df$Sample_ID, substring(rownames(table_cluster_assignment), 1, 4))]
df$ICRscore = table_cluster_assignment$ICRscore[match(df$Sample_ID, substring(rownames(table_cluster_assignment), 1, 4))]
df$histology = clinical_data$Tumor_morphology[match(substring(df$Sample_ID, 1, 3), clinical_data$Patient_ID)]

df$Mutation_cat = NA
df$Mutation_cat[which(df$Mutational_load_per_MB <= 12)] = "nonhypermutated"
df$Mutation_cat[which(df$Mutational_load_per_MB > 12)] = "hypermutated"

df$ICR = factor(df$ICR, levels = c("ICR Low", "ICR Medium", "ICR High"))
df$CMS = factor(df$CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
df$MSI = factor(df$MSI, levels = c("MSI-H", "MSS"))
df$Mutation_cat = factor(df$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))

# Define order based on mutational load
df = df[order(df$ICRscore),]
df = df[order(df$Mutational_load),]
df = df[order(df$CMS),]

if(version %in% c("G", "H")){
  df = df[order(df$Mutation_cat),]
  df = df[order(df$ICR),]
}else{
  df = df[order(df$ICR),]
  df = df[order(df$Mutation_cat),]
}


if(subset == "all_patients"){}else{
  df = df[which(df$Mutation_cat == subset),]
}
sample_order = rownames(df)
mat = mat[,sample_order]

# Heatmap annotation
ha = HeatmapAnnotation(`Mutational load` = anno_barplot(df$Mutational_load, baseline = 0,
                                                        axis_param = list(
                                                          side = "right"
                                                        )),
                       `Mutation status` = df$Mutation_cat,
                       ICR = df$ICR,
                       CMS = df$CMS,
                       MSI = df$MSI,
                       Histology = df$histology,
                       annotation_name_side = "left",
                       col = list(`Mutation status` = c("nonhypermutated" = "#A7EABD", "hypermutated" = "#EAAED0"),
                                  ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue"),
                                  CMS = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
                                          "mixed" = "lightgrey"),
                                  MSI = c("MSI-H" = "purple", "MSS" = "yellow"),
                                  Histology = c("adenocarcinoma nno" = "grey", "adenocarcinoom in villeus adenoom" = "green", 
                                                "blue", "adenocarcinoom, intestinaal type" = "purple", "cribriform carcinoom" = "red", 
                                                "mucineus adenocarcinoom" = "darkorange",
                                                "zegelringcel carcinoom" = "darkgreen", "cyan", "adenocarcinoom met gemengde subtypes" = "darkblue", "neoplasma, maligne" = "pink")
                       )
)

# Plotting
dir.create("./Figures/WES/009_Oncoprint", showWarnings = FALSE)
pdf(paste0("./Figures/WES/009_Oncoprint/Fig3d_009_", version, "_OncoPrint.pdf"), width = 11, height = 7) # 11 #7.5,  #20 #18
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          top_annotation = ha,
          row_names_max_width = unit(12, "in"),
          row_order = rownames(mat),
          column_order = sample_order)
dev.off()

