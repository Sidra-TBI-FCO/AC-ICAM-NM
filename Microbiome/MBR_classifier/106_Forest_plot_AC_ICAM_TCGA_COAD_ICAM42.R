
# Validation of MBR from AC-ICAM in ICAM42, TCGA-COAD, and ICAM42+TCGA-COAD

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx", "stringr", "forestplot"))

# set parameters
survival = "OS" # OS

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/105_Collect_Data_for_forest_plot/105_HR_table_AC_ICAM_ICAM42_TCGA_COAD.Rdata")
results_df = HR_table_all

#remove digits using (signif)
results_df$p_val = signif(results_df$p_val, digits=3)
results_df$HR = signif(results_df$HR, digits=3)
results_df$p_val_logrank = signif(results_df$p_val_logrank, digits=3)

results_df$Survival = gsub("TCGA_COAD", "TCGA-COAD", results_df$Survival)
results_df$Survival = gsub("ICAM42_", "ICAM42 & ", results_df$Survival)
  

write.csv(results_df, file ="./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/105_Collect_Data_for_forest_plot/105_HR_table_AC_ICAM_ICAM42_TCGA_COAD.csv",
          row.names = FALSE)
# Order based on HR column
#results_df = results_df[order(results_df$HR),]

# set t.test_results$pathway as length 
N.sig = length(results_df$Survival)

#x = N.sig + 2
x = N.sig

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,results_df$HR[1:x])),
    lower = c(NA,results_df$CI_lower[c(1:x)]),
    upper = c(NA,results_df$CI_upper[c(1:x)])),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

tabletext<-cbind(
  c("Survival", as.character(results_df$Survival)[c(1:x)]),
  #c("p-value-logrank", results_df$p_val_logrank[c(1:x)]),
  c("p-value", results_df$p_val[c(1:x)]),
  c("HR",      results_df$HR[c(1:x)]))

svg(filename = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/forestplot/foresplot_AC-ICAM_TCGA_COAD_ICAM42.svg"), 
    height = 2, width = 5)
forestplot(tabletext, 
           graph.pos = 2, # position of the graph 2 (after first column)
           cochrane_from_rmeta,new_page = FALSE,
           is.summary=c(TRUE,rep(FALSE,x),TRUE,rep(FALSE,x),TRUE,FALSE),
           #clip=c(0,8),
           xlog=T, # If TRUE, middle line will be at 1, if false middle line will be at 0 
           xticks = c(0.05,0.2,1),
           boxsize = .15, #.25
           vertices = TRUE, # Set this to TRUE if you want the ends of the confidence intervals to be shaped as a T
           col=fpColors(box="black",line="black", summary="black"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10), cex = 1, ticks = gpar(fontsize = 20))) #13
dev.off()

#######

pdf(file = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/forestplot/forestplot_AC-ICAM_TCGA_COAD_ICAM42.pdf"),height = 3, width = 5)
forestplot(tabletext, 
           graph.pos = 2, # position of the graph 2 (after first column)
           cochrane_from_rmeta,new_page = FALSE,
           is.summary=c(TRUE,rep(FALSE,x),TRUE,rep(FALSE,x),TRUE,FALSE),
           #clip=c(0,8),
           xlog=T, # If TRUE, middle line will be at 1, if false middle line will be at 0 
           xticks = c(0.05,0.2,1),
           boxsize = .15, #.25
           vertices = TRUE, # Set this to TRUE if you want the ends of the confidence intervals to be shaped as a T
           col=fpColors(box="black",line="black", summary="black"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10), cex = 1, ticks = gpar(fontsize = 20))) #13 
dev.off()

