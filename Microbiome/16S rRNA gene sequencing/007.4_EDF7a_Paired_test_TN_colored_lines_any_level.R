
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix", "scales"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"


# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))

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
                #ICR_cluster = NA,
                stringsAsFactors = FALSE)

dir.create("./Figures/Microbiome/007.4_Paired_wilcoxon_test_TN", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/007.4_Paired_wilcoxon_test_TN/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/007.4_Paired_wilcoxon_test_TN/", Type, "/", Rank), showWarnings = FALSE)

i = 134 #134 #419 #177
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$T = abundance_T[,micr][match(df$Patient_ID, rownames(abundance_T))]
  df$N = abundance_N[,micr][match(df$Patient_ID, rownames(abundance_N))]
  df$Delta = df$T - df$N
  
  if(Type == "Relative"){
    df$T[which(df$T == 0)] = 1e-4
    df$N[which(df$N == 0)] = 1e-4
  }
  
  melt = melt(df)
  melt$variable = factor(melt$variable, levels = c("N", "T"))
  #levels(melt$variable) = c("Normal", "Tumor")
  mean = mean(melt$value)
  melt$Delta = df$Delta[match(melt$Patient_ID, df$Patient_ID)]
  melt = melt[-which(is.na(melt$variable)),]
  
  melt$Delta2 = NA
  melt$Delta2[which(melt$Delta < 0)] = "down"
  melt$Delta2[which(melt$Delta == 0)] = "same"
  melt$Delta2[which(melt$Delta > 0)] = "up"
  
  ylabel = gsub(".*\\D_5__", "", micr)
  
  if(make_plot == "make_plot" | mean > 0.001){
    plot = ggplot(melt, aes(x = variable, y = value)) +
      #facet_grid(cols = vars(ICR_cluster)) +
      geom_boxplot(outlier.shape = NA) +
      #geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
      geom_line(aes(group = Patient_ID, color = Delta2)) +
      scale_colour_manual(values = c("down" = alpha("blue", 0.2), "same" = alpha("black", 0.2), "up" = alpha("red", 0.2))) +
      geom_point(size = 0.8) +
      ylab(paste0(Type, " abundance \n", ylabel)) +
      theme_bw() +
      scale_y_log10(labels = comma) +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black"),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            strip.background = element_blank(),
            strip.text = element_text(size = 15, colour = "black"),
            legend.position = "none")
    
    pdf(paste0("./Figures/Microbiome/007.4_Paired_wilcoxon_test_TN/", Type, "/", Rank, "/EDF7a_",
               micr, "_paired_boxplot_N_T.pdf"), width = 3, height = 3.5)
    plot(plot)
    dev.off()
  }
}
