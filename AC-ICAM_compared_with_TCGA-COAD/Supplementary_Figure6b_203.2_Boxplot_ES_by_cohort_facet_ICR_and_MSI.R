

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "easyGgplot2"))

dir.create("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC", showWarnings = FALSE)
dir.create("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003.2_Boxplot_ES_by_cohort", showWarnings = FALSE)

# Set parameters
Geneset = "ConsensusTME" # "ConsensusTME" or "Hallmarks"
Signatures = c("B_cells", "T_cells_CD8", "NK_cells", "ICRscore")
subgroup = "ICR_cluster" # "MSI" or "ICR_cluster"

# Load data
if(Geneset == "ConsensusTME"){
  load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ssGSEA/ConsensusTME_TCGA_COAD_Sidra_LUMC.Rdata")
  load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ICRscores_TCGA_AC_ICAM.Rdata")
  ES = rbind(ES, t(ICRscores))
  ES = ES[Signatures,]
}
load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ICR_clusters_rbind.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
MANTIS_AC_ICAM = MANTIS
load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/External_Data/From_TCGA_Masterfile/TCGA_MANTIS.Rdata")
MANTIS_TCGA_COAD = MANTIS

MANTIS_AC_ICAM = MANTIS_AC_ICAM[, colnames(MANTIS_TCGA_COAD)]
MANTIS = rbind(MANTIS_AC_ICAM, MANTIS_TCGA_COAD)

# Prepare data for plot
df = data.frame(Sample_ID = colnames(ES), ES = NA, cohort = "AC-ICAM")

df$cohort[grep("TCGA", df$Sample_ID)] = "TCGA-COAD"
df$ICR_cluster = table_cluster_assignment$ICR_HML[match(df$Sample_ID, rownames(table_cluster_assignment))]
table(df$ICR_cluster, exclude = NULL)
df$ICR_cluster = factor(df$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
df$Patient_ID = substring(df$Sample_ID, 1, 12)
df$Patient_ID[which(df$cohort == "AC-ICAM")] = substring(df$Sample_ID[which(df$cohort == "AC-ICAM")], 1, 3)
df$MSI = MANTIS$MSI[match(df$Patient_ID, MANTIS$Patient_ID)]
table(df$MSI, exclude = NULL)
df$MSI = factor(df$MSI, levels = c("MSS", "MSI-H"))

if(subgroup == "ICR_cluster"){
  df = df[-which(df$ICR_cluster == "ICR Medium"),]
}

if(subgroup == "MSI"){
  df = df[-which(is.na(df$MSI)),]
}

cell_subsets = rownames(ES)

results = data.frame(Signature = cell_subsets, p_value = NA, FDR = NA, mean_AC_ICAM = NA, mean_TCGA = NA)

i=2
for (i in 1:length(cell_subsets)){
  subset = cell_subsets[i]
  df$ES = ES[subset,][match(df$Sample_ID, colnames(ES))]
  
  plot = ggplot(df, aes(x = cohort, y = ES)) +
    facet_grid(.~get(subgroup)) +
    geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
    scale_fill_manual(values = c("AC-ICAM" = "#E7A686", "TCGA-COAD" = "#C5EBED")) +
    geom_jitter(width = 0.2, size = 0.8) +
    stat_compare_means(method = "t.test", 
                       comparisons = list(c("AC-ICAM", "TCGA-COAD"))) +
    ylab(paste0(gsub("_", " ", subset), "\nenrichment score")) +
    xlab("") +
    theme_bw() +
    theme(axis.title.x = element_text(colour = "black", size = 15),
          axis.title.y = element_text(colour = "black", size = 15),
          axis.text.x = element_text(colour = "black", size = 15, angle = 45,
                                     vjust = 1, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 15),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(colour = "black", size = 15))
  
  assign(paste0("p", i, sep = ""), plot)
  
  png(paste0("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003.2_Boxplot_ES_by_cohort/",Geneset, "_", subset, "_ES_by_", subgroup,
             "_by_cohort.png"),
      res = 600, units = "in", width = 4, height = 5)
  plot(plot)
  dev.off()
  
  #df = df[which(df$ICR_cluster == "ICR Low"),]
  #test = t.test(df$ES[which(df$cohort == "AC-ICAM")], df$ES[which(df$cohort == "TCGA-COAD")])
  #test$p.value
  #results$p_value[which(results$Signature == subset)] = test$p.value
  #results$mean_AC_ICAM[which(results$Signature == subset)] = test$estimate[1]
  #results$mean_TCGA[which(results$Signature == subset)] = test$estimate[2]
}

plots = paste("p", 1:length(cell_subsets), sep = "")
list_of_plots = mget(plots)

results$FDR = p.adjust(results$p_value, method = "BH", n = nrow(results))
results$Delta = results$mean_TCGA - results$mean_AC_ICAM

if(Geneset == "ConsensusTME" ){
  pdf(paste0("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003.2_Boxplot_ES_by_cohort/Supplementary_Figure6b_", subgroup, "_", Geneset,
             "_", subgroup,".pdf"),
      width = 15, height = 6)
  ggplot2.multiplot(plotlist = list_of_plots[1:length(cell_subsets)], cols=4)
  dev.off()
}




