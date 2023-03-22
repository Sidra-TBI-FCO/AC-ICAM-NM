
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr", "stringr"))

# Set parameters


# Load data
load("./Processed_Data/Microbiome/23_Alpha_diversity/Sample_based_filted/alpha_diversity_TN_pairs.Rdata") 

dir.create("./Figures/Microbiome/024.2_Boxplot_alpha_diversity_by_TN", showWarnings = FALSE)

diversity_indices = colnames(df_alpha[4:7])

test_df = data.frame(diversity_index = diversity_indices, p_val = NA)

i=1
for (i in 1:length(diversity_indices)){
  index = diversity_indices[i]
  
  df_plot = df_alpha[, c("Patient_ID", "Tissue", index)]
  
  plot = ggplot(df_plot, aes(x = Tissue, y = get(index))) +
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
    geom_point(size = 0.8) +
    theme_bw() +
    xlab("") +
    ylab(index) +
    #stat_compare_means(method = "wilcox.test", comparisons =  list(c("N", "T"))) +
    theme(axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 15, color = "black"),
          aspect.ratio = 1.8/1)
  
  wil = wilcox.test(df_plot[, index][which(df_plot$Tissue == "T")], 
                    df_plot[, index][which(df_plot$Tissue == "N")],
                    paired = TRUE)
 
  test_df$p_val[which(test_df$diversity_index == index)] = wil$p.value
  
  png(paste0("./Figures/Microbiome/024.2_Boxplot_alpha_diversity_by_TN/Boxplot_", index, "_by_TN.png"),
      res = 600, width = 2.3, height = 3.3, units = "in")
  plot(plot)
  dev.off()
  
  pdf(paste0("./Figures/Microbiome/024.2_Boxplot_alpha_diversity_by_TN/Boxplot_", index, "_by_TN.pdf"),
      width = 2.3, height = 3.3)
  plot(plot)
  dev.off()
}

