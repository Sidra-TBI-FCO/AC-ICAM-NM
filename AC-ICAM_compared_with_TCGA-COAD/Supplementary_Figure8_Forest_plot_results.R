
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

# Set parameters
Random = "Random" # "Random" or ""

# Load data
load(paste0("./Analysis/Exploration_reviewers/Purity_cutoff_ICR_Survival/Results_", Random, 
            "_Survival_Analysis_ICRscore_continuous.Rdata"))

HR_table = results_df

## Forest plot seperate script
n_signatures = nrow(HR_table)
x = n_signatures + 2

HR_table = HR_table[order(HR_table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_signatures]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_signatures)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_signatures)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

HR_table$p_value = signif(HR_table$p_value, 2)
HR_table$HR = signif(HR_table$HR, 3)
tabletext<-cbind(
  c("p-value", HR_table$p_value[c(1:n_signatures)]),
  c("HR",      HR_table$HR[c(1:n_signatures)]))

#dev.new()
dir.create("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Forest_plots", showWarnings = FALSE)

pdf(file = paste0("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Forest_plots/Sup_Figure8_Forest_plot_",
                  Random, "_ICRscore_OS_Survival.pdf"),
    height = 25, width = 15)
forestplot(tabletext,
           cochrane_from_rmeta,new_page = FALSE,
           is.summary=c(TRUE,rep(FALSE,n_signatures),TRUE,rep(FALSE,n_signatures),TRUE,FALSE),
           #clip=c(0,8),
           xlog=TRUE,
           xlim = c(-2, 2),
           boxsize = .25,
           cex = 1.5,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 10,cex=1.5), cex = 1.5))
dev.off()
