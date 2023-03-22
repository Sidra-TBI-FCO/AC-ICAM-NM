
### Association number of templates and ICR score

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
data_type = "raw" # "raw" or "DNA_input_corrected"

# Set parameters and load data
if(data_type == "raw"){
  load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
  variables = colnames(TCR_Overview)[c(2:8, 23)]
}
if(data_type == "DNA_input_corrected"){
  load("./Analysis/TCR/2_Corrections_for_DNA_input/TCR_Overview_DNA_input_corrected.Rdata")
  variables = c("T_cell_fraction_in_percent", "productive_templates_cor", "total_rearrangements_cor", "productive_rearrangements_cor")
}

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]

N.variables = length(variables)

i = 6
for (i in 1:N.variables){
  variable = variables[i]
  # Append ICR Cluster data to the overview file
 
  TCR_Overview[,variable] = as.numeric(TCR_Overview[,variable])
  
  plot = ggplot(TCR_Overview, aes(x = ICRscore, y = TCR_Overview[,variable])) +
    geom_point(aes(color = HLM_cluster), size = 0.8) +
    ylab(gsub("\\_", " ", variable)) +
    scale_color_manual(values = c("red", "green", "blue")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          aspect.ratio = 1/1) +
    labs(color = "ICR cluster") +
    #stat_cor(method = "pearson", size = 5) +
    geom_smooth(method="lm")
    
  
  dir.create("./Figures", showWarnings = FALSE)
  dir.create("./Figures/TCR", showWarnings = FALSE)
  dir.create("./Figures/TCR/1_Scatterplots_ICR_TCRmetrics", showWarnings = FALSE)
  pdf(paste0("./Figures/TCR/1_Scatterplots_ICR_TCRmetrics/1_v3_Scatterplot_ICR_", variable, ".pdf"),
      width = 5, height = 3)
  plot(plot)
  dev.off()
}
