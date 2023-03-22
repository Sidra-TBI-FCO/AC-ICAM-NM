
### T cell metrics

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Load data
load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

ES = t(immune_sig_df)
ES = as.data.frame(ES)

ES = as.matrix(ES)
TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]

variables = rownames(ES)
N.variables = length(variables)

results_ICR = data.frame(Variable = variables, cor_coef = NA)
results_productive_templates = data.frame(Variable = variables, cor_coef = NA, p_val = NA)
results_productive_templates_cor = data.frame(Variable = variables, cor_coef = NA)
results_clonality = data.frame(Variable = variables, cor_coef = NA, p_val = NA)

i = 1
for (i in 1:N.variables){
  variable = variables[i]
  # Append ICR Cluster data to the overview file
  TCR_Overview[,variable] = ES[variable, ][match(TCR_Overview$sample_name, substring(colnames(ES), 1, 4))]
  TCR_Overview[,variable] = as.numeric(TCR_Overview[,variable])
  
  correlation = cor(x = TCR_Overview$ICRscore, y = TCR_Overview[,variable], method = "pearson")
  results_ICR$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  plot2 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_templates)) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    xlab(gsub("\\_", " ", variable)) +
    ylab("productive templates") +
    scale_color_manual(values = c("red", "green", "blue")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          legend.position = "none",
          aspect.ratio = 1/1) +
    labs(color = "ICR cluster") +
    stat_cor(method = "pearson", size = 5) +
    geom_smooth(method="lm")
  
  cor_test = cor.test(x = TCR_Overview[,variable], y = TCR_Overview$productive_templates, method = "pearson")
  correlation = cor_test$estimate
  p_val = cor_test$p.value
  results_productive_templates$cor_coef[which(results_ICR$Variable == variable)] = correlation
  results_productive_templates$p_val[which(results_ICR$Variable == variable)] = p_val
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  #png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_v3_Scatterplot_", variable, "_TCR_productive_templates.png"), res = 600,
   #   width = 3.2, height = 3, units = "in")
  #plot(plot2)
  #dev.off()
  
  correlation = cor(x = TCR_Overview[,variable], y = TCR_Overview$productive_templates_cor, method = "pearson")
  results_productive_templates_cor$cor_coef[which(results_ICR$Variable == variable)] = correlation
  
  plot4 = ggplot(TCR_Overview, aes(x = TCR_Overview[,variable], y = productive_clonality)) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    xlab(gsub("\\_", " ", variable)) +
    ylab("productive clonality") +
    scale_color_manual(values = c("red", "green", "blue")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          legend.position = "none",
          aspect.ratio = 1/1) +
    labs(color = "ICR cluster") +
    stat_cor(method = "pearson", size = 5) +
    geom_smooth(method="lm") 
  #xlab("") +
  #ylab("")
  
  cor_test = cor.test(x = TCR_Overview[,variable], y = TCR_Overview$productive_clonality, method = "pearson")
  correlation = cor_test$estimate
  p_val = cor_test$p.value
  results_clonality$cor_coef[which(results_ICR$Variable == variable)] = correlation
  results_clonality$p_val[which(results_ICR$Variable == variable)] = p_val
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR/6_Scatterplots_T_cell_metrics", showWarnings = FALSE)
  #png(paste0("./Figures/TCR/6_Scatterplots_T_cell_metrics/6_v3_Scatterplot_", variable, "_TCR_cell_clonality.png"), res = 600,
   #   width = 3.2, height = 3, units = "in")
  #plot(plot4)
  #dev.off()
}

results_ICR = results_ICR[order(results_ICR$cor_coef, decreasing = TRUE),]
results_productive_templates = results_productive_templates[order(results_productive_templates$cor_coef, decreasing = TRUE),]
results_clonality = results_clonality[order(results_clonality$cor_coef, decreasing = TRUE),]

results_ICR$cor_coef = round(results_ICR$cor_coef, 3)
results_productive_templates$cor_coef = round(results_productive_templates$cor_coef, 3)
results_clonality$cor_coef = round(results_clonality$cor_coef, 3)

results_productive_templates$p_val = signif(results_productive_templates$p_val, 3)
results_clonality$p_val = signif(results_clonality$p_val, 3)

results_productive_templates$FDR = p.adjust(results_productive_templates$p_val, method = "BH")
results_clonality$FDR = p.adjust(results_clonality$p_val, method = "BH")

dir.create("./Analysis/TCR/006_T_cell_infiltration_cor", showWarnings = FALSE)
write.csv(results_ICR, file = "./Analysis/TCR/006_T_cell_infiltration_cor/April_2021_All_clean_results_ICR.csv", row.names = FALSE)
write.csv(results_productive_templates, file = "./Analysis/TCR/006_T_cell_infiltration_cor/April_2021_All_clean_results_productive_templates.csv", row.names = FALSE)
#write.csv(results_productive_templates_cor, file = paste0("./Analysis/TCR/006_All_T_cell_infiltration_cor/", Source, "results_productive_templates_cor.csv"), row.names = FALSE)
write.csv(results_clonality, file = "./Analysis/TCR/006_T_cell_infiltration_cor/April_2021_All_clean_results_clonality.csv", row.names = FALSE)
