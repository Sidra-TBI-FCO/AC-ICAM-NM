# common samples OS

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx", "stringr", "forestplot"))

# set parameters
survival = "DFS" # OS  DFS 
matrix = "RA" # RA # OTU
version = "without" # without # subtract

# load data
load(paste0("./Analysis/Exploration_reviewers/Microbiome/round_two/SparCC/OS/sparcc_HR_table_",survival, "_network_clusters_",version,"_Prevotella2_",matrix,".Rdata"))


#remove digits using (signif)
results_df = results
results_df$cox_pval = signif(results_df$cox_pval, digits=2)
results_df$HR = signif(results_df$HR, digits=2)
results_df$logrank_pval = signif(results_df$logrank_pval, digits=2)

results_df$FDR = p.adjust(results_df$cox_pval, method = "fdr", n = nrow(results_df))
results_df$FDR = signif(results_df$FDR, digits=2)

results_df$Cluster = gsub("Group", "Group ", results_df$Cluster)

# Order based on HR column
#results_df = results_df[order(results_df$HR),]

# set t.test_results$pathway as length 
N.sig = length(results_df$Cluster)

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
  c("Survival", as.character(results_df$Cluster)[c(1:x)]),
  #c("p-value-logrank", results_df$p_val_logrank[c(1:x)]),
  c("p-value", results_df$cox_pval[c(1:x)]),
  c("FDR", results_df$FDR[c(1:x)]),
  c("HR",  results_df$HR[c(1:x)]))

# svg(filename = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/forestplot/foresplot_",survival,"_",cohort,"_samples_all_platforms_updated.svg"), height = 3, width = 5)
# forestplot(tabletext, 
#            graph.pos = 2, # position of the graph 2 (after first column)
#            cochrane_from_rmeta,new_page = FALSE,
#            is.summary=c(TRUE,rep(FALSE,x),TRUE,rep(FALSE,x),TRUE,FALSE),
#            #clip=c(0,8),
#            xlog=T, # If TRUE, middle line will be at 1, if false middle line will be at 0 
#            xticks = c(0.05,0.2,1),
#            boxsize = .15, #.25
#            vertices = TRUE, # Set this to TRUE if you want the ends of the confidence intervals to be shaped as a T
#            col=fpColors(box="black",line="black", summary="black"),
#            txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10), cex = 1, ticks = gpar(fontsize = 20))) #13
# dev.off()

#######

pdf(file = paste0("./Figures/Exploration_reviewers/Microbiome/round_2/SparCC/OS/",survival,"_forestplot_network_clusters_",version,"_Prevotella2_",matrix,".pdf"),height = 2, width = 5)

forestplot(tabletext,
           graph.pos = 2, # position of the graph 2 (after first column)
           cochrane_from_rmeta,new_page = FALSE,
           is.summary=c(TRUE,rep(FALSE,x),TRUE,rep(FALSE,x),TRUE,FALSE),
           #clip=c(0,8),
           xlog=T, # If TRUE, middle line will be at 1, if false middle line will be at 0
           xticks = c(0.2, 0.5,1, 1.5, 2),
           boxsize = .15, #.25
           vertices = TRUE, # Set this to TRUE if you want the ends of the confidence intervals to be shaped as a T
           col=fpColors(box="black",line="black", summary="black"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10), cex = 1, ticks = gpar(fontsize = 20))) #13
dev.off()
