

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.3_147_Genera_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_or_N_246_samples.Rdata")

rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                T = NA,
                N = NA,
                ratio = NA,
                delta = NA,
                stringsAsFactors = FALSE)

results = data.frame(Name = rownames(rank_abundancies), p_val = NA, wilcox_statistic = NA,
                     mean_N = NA, mean_T = NA, median_N = NA, median_T = NA, 
                     mean_ratio_paired = NA, mean_delta_paired = NA)

dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN", showWarnings = FALSE)
dir.create(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/", Type, "/", Rank), showWarnings = FALSE)

max_new = 1
min_new = 1

i = 14
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$T = abundance_T[,micr][match(df$Patient_ID, rownames(abundance_T))]
  df$N = abundance_N[,micr][match(df$Patient_ID, rownames(abundance_N))]
  
  df2 = df
  df2$T_man_cor = df2$T + 1e-4   # manual correction for abundancies
  df2$N_man_cor = df2$N + 1e-4
  df2$ratio_man_cor = df2$T_man_cor / df2$N_man_cor # check
  df2$ratio_man_cor[which(df2$T == df2$N)] = 1 # check
  
  df2$ratio = df2$T / df2$N
  df2$ratio[which(df2$T == df2$N)] = 1 # check
  df2$ratio[which(df2$ratio == Inf & df2$T > df2$N)] = 30  # value of max_new maximum fold change assigned to the ones coming from zero
  df2$ratio[which(df2$ratio == 0)] = 0.01  # value of min_new maximum fold change assigned to the ones coming from zero
  df2$delta = df2$T - df2$N
  max = max(df2$ratio)
  min = min(df2$ratio)
  if(max > max_new){
    max_new = max
  }
  if(min < min_new){
    min_new = min
  }
  
  test = wilcox.test(x = df$N, y = df$T, paired = TRUE, alternative = "two.sided")
  results$p_val[which(results$Name == micr)] = test$p.value
  results$wilcox_statistic[which(results$Name == micr)] = test$statistic
  results$mean_T[which(results$Name == micr)] = mean(df$T)
  results$mean_N[which(results$Name == micr)] = mean(df$N)
  results$median_T[which(results$Name == micr)] = median(df$T)
  results$median_N[which(results$Name == micr)] = median(df$N)
  results$mean_ratio_paired[which(results$Name == micr)] = mean(df$ratio, na.rm = TRUE) # check
  results$mean_delta_paired[which(results$Name == micr)] = mean(df$delta, na.rm = TRUE)
  
  if(Type == "Relative"){
    df$T[which(df$T == 0)] = 1e-4
    df$N[which(df$N == 0)] = 1e-4
  }
  
  melt = melt(df)
  melt$variable = factor(melt$variable, levels = c("N", "T"))
  levels(melt$variable) = c("Normal", "Tumor")
  
  if(test$statistic == 0){next}
  if(test$p.value < 0.05 | make_plot == "make_plot"){
    plot = ggplot(melt, aes(x = variable, y = value)) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
      geom_point(size = 0.8) +
      scale_y_log10() +
      ylab(paste0(micr, "\n", Type, " abundance")) +
      theme_bw() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black"),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"))
    
    #png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/", Type, "/", Rank, "/",
    #          micr, "_paired_boxplot_N_T.png"), res = 600, units = "in", width = 3, height = 3.5)
    #plot(plot)
    #dev.off()
  }
}

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))
results = results[order(results$p_val),]

results$Direction = NA
results$Direction[which(results$mean_T > results$mean_N)] = "Enriched in tumor"
results$Direction[which(results$mean_T <= results$mean_N)] = "Enriched in normal"
results$Ratio_mean = results$mean_T / results$mean_N

dir.create(paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_Wilcoxon_test"), showWarnings = FALSE)

save(results, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_Wilcoxon_test/48_TN_Paired_Wilcoxon_test.Rdata"))
write.csv(results, file = paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_Wilcoxon_test/48_TN_Paired_Wilcoxon_test.csv"),
          row.names = FALSE)
