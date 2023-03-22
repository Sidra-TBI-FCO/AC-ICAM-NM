
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "easyGgplot2"))

dir.create("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC", showWarnings = FALSE)
dir.create("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003_Boxplot_ES_by_cohort", showWarnings = FALSE)

# Set parameters
Geneset = "Hallmarks" # "ConsensusTME" or "Hallmarks"

# Load data
if(Geneset == "ConsensusTME"){
  load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ssGSEA/ConsensusTME_TCGA_COAD_Sidra_LUMC.Rdata")
  load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ICRscores_TCGA_AC_ICAM.Rdata")
  ES = rbind(ES, t(ICRscores))
}
if(Geneset == "Hallmarks"){
  load("./Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ssGSEA/Hallmark_pathways_TCGA_COAD_Sidra_LUMC.Rdata")
}

# Prepare data for plot
df = data.frame(Sample_ID = colnames(ES), ES = NA, cohort = "AC-ICAM")

df$cohort[grep("TCGA", df$Sample_ID)] = "TCGA-COAD"

cell_subsets = rownames(ES)

results = data.frame(Signature = cell_subsets, p_value = NA, FDR = NA, mean_AC_ICAM = NA, mean_TCGA = NA)

i=1
for (i in 1:length(cell_subsets)){
  subset = cell_subsets[i]
  df$ES = ES[subset,][match(df$Sample_ID, colnames(ES))]
  
  plot = ggplot(df, aes(x = cohort, y = ES)) +
    geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
    scale_fill_manual(values = c("AC-ICAM" = "#E7A686", "TCGA-COAD" = "#C5EBED")) +
    geom_jitter(width = 0.2, size = 0.8) +
    stat_compare_means(method = "t.test", label = "p.signif",
                       comparisons = list(c("AC-ICAM", "TCGA-COAD"))) +
    ylab(paste0(gsub("_", " ", subset), "\nenrichment score")) +
    xlab("") +
    theme_bw() +
    theme(axis.title.x = element_text(colour = "black", size = 15),
          axis.title.y = element_text(colour = "black", size = 15),
          axis.text.x = element_text(colour = "black", size = 15),
          axis.text.y = element_text(colour = "black", size = 15),
          legend.position = "none")
  
  assign(paste0("p", i, sep = ""), plot)
  
  png(paste0("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003_Boxplot_ES_by_cohort/",Geneset, "_", subset, "_ES_by_cohort.png"),
      res = 600, units = "in", width = 4, height = 5)
  plot(plot)
  dev.off()
  
  test = t.test(df$ES[which(df$cohort == "AC-ICAM")], df$ES[which(df$cohort == "TCGA-COAD")])
  results$p_value[which(results$Signature == subset)] = test$p.value
  results$mean_AC_ICAM[which(results$Signature == subset)] = test$estimate[1]
  results$mean_TCGA[which(results$Signature == subset)] = test$estimate[2]
}

plots = paste("p", 1:length(cell_subsets), sep = "")
list_of_plots = mget(plots)

results$FDR = p.adjust(results$p_value, method = "BH", n = nrow(results))
results$Delta = results$mean_TCGA - results$mean_AC_ICAM

if(Geneset == "ConsensusTME" ){
  pdf(paste0("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003_Boxplot_ES_by_cohort/", Geneset,".pdf"),
      width = 15, height = 20)
  ggplot2.multiplot(plotlist = list_of_plots[1:length(cell_subsets)], cols=4)
  dev.off()
}

if(Geneset == "Hallmarks"){
  results = results[-which(results$Signature == "HALLMARK_SPERMATOGENESIS"),]
  results$FDR = p.adjust(results$p_value, method = "BH", n = nrow(results))
  index = which(results$Delta > 0 & results$p_value < 0.01)
  list_of_plots = list_of_plots[index]
  
  pdf(paste0("./Figures/030_Merged_TCGA_COAD_Sidra_LUMC/003_Boxplot_ES_by_cohort/", Geneset,".pdf"),
      width = 15, height = 13)
  ggplot2.multiplot(plotlist = list_of_plots[1:length(list_of_plots)], cols=4)
  dev.off()
}



dir.create("./Analysis/Exploration_reviewers/ES_by_cohort", showWarnings = FALSE)
write.csv(results, file = paste0("./Analysis/Exploration_reviewers/ES_by_cohort/Supplementary_Figure5_ES_", Geneset, "_by_cohort.csv"),
          row.names = FALSE)



