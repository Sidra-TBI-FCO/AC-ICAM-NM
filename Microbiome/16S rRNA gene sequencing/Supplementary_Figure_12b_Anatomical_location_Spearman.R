
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "plotrix", "RColorBrewer"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
Tissue = "N" # "T" or "N" or "Difference"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
make_plot = "do_not_make_plot" # "do_not_make_plot" "make_plot"

# Load data
if(Tissue == "T"){
  load(paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_", 
              "T_246_samples.Rdata"))
}
if(Tissue == "N"){
  load(paste0("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_N_246_samples_based_on_normal.Rdata"))
}

load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)

abundance = abundance_T

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                Abundance = NA,
                Anatomical_location = NA,
                Primary_tumor_side = NA,
                stringsAsFactors = FALSE)

df$Anatomical_location = clinical_data$tumour_anatomic_site[match(df$Patient_ID, clinical_data$Patient_ID)]
df$Primary_tumor_side[which(df$Anatomical_location %in% c("ceceum", "colon ascendens", "flexura hepatica", 
                                                          "colon transversum"))] = "Right sided"
df$Primary_tumor_side[which(df$Anatomical_location %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                          "rectosigmoideum"))] = "Left sided"
df$Primary_tumor_side = factor(df$Primary_tumor_side, levels = c("Right sided", "Left sided"))

df$Anatomical_location[which(df$Anatomical_location == "ceceum")] = "caecum"
df$Anatomical_location = factor(df$Anatomical_location, levels = c("caecum", "colon ascendens", "flexura hepatica", 
                                                                   "colon transversum", "flexura lienalis", "colon descendens", 
                                                                   "colon sigmoideum", "rectosigmoideum"))

table(df$Anatomical_location)
levels(df$Anatomical_location) = c("1", "2", "3", "4", "5", "6", "7", "8")
df$Anatomical_location = as.numeric(df$Anatomical_location)
table(df$Anatomical_location)

results = data.frame(Name = rownames(rank_abundancies), cor = NA, p_val = NA)

dir.create("./Figures/Microbiome/Re_buttal_014.2_Spearman", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/", Type), showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/", Type, "/", Rank), showWarnings = FALSE)

i = 105  # 137/128 "Akkermansia" or 93/90 "NK4A214" or 67/67 "UCG-004" or 108/105 "UCG-003"
for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[,micr][match(df$Patient_ID, rownames(abundance))]
  
  cor = cor.test(df$Abundance, df$Anatomical_location, method = "spearman")
  results$p_val[which(results$Name == micr)] = cor$p.value
  results$cor[which(results$Name == micr)] = cor$estimate
  
  if(Type == "Relative" & Tissue %in% c("T", "N")){
    df$Abundance[which(df$Abundance == 0)] = 1e-4
  }
  
  melt = df
  
  if(is.na(cor$p.value)){next}
  if(cor$p.value < 0.05 | make_plot == "make_plot"){
    melt$Anatomical_location = factor(melt$Anatomical_location)
    levels(melt$Anatomical_location) = c("ceceum", "colon ascendens", "flexura hepatica", 
                                         "colon transversum", "flexura lienalis", "colon descendens", 
                                         "colon sigmoideum", "rectosigmoideum")
    label = gsub(".*\\D_5__", "", micr)
    
    plot = ggplot(melt, aes(x = Anatomical_location, y = Abundance, fill = Anatomical_location)) +
      scale_fill_manual(values = brewer.pal(n = 8, name = "Set3")) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(aes(group = Patient_ID), color = alpha("lightgrey", 0.5)) +
      geom_point(size = 0.8) + 
      ylab(paste0(label)) +
      theme_bw() +
      scale_y_log10() +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(size = 15, colour = "black",
                                       angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            legend.position = "none")
    
    pdf(paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/", Type, "/", Rank, 
               "/Supplementary_Figure_12b_", Tissue, "_",
               micr, "_boxplot_anatomic location.pdf"), width = 3, height = 3.5)
    plot(plot)
    dev.off()
  }
}

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))

results = results[order(results$p_val),]

dir.create("./Analysis/Microbiome/Re_buttal_014.2_Spearman", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/Re_buttal_014.2_Spearman/",  Type), showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/Re_buttal_014.2_Spearman/",  Type, "/", Rank), showWarnings = FALSE)

#save(results, file = paste0("./Analysis/Microbiome/Re_buttal_014.2_Spearman/",  Type, "/", Rank,
#                            "/Spearman_", Type, "_", Rank, "_", Tissue, ".Rdata"))
#write.csv(results, file = paste0("./Analysis/Microbiome/Re_buttal_014.2_Spearman/",  Type, "/", Rank,
 #                                "/Spearman_", Type, "_", Rank, "_", Tissue, ".csv"),
  #        row.names = FALSE)

results_FDR = results[which(results$FDR < 0.1),]
results_FDR$direction = "Lower"
results_FDR$direction[which(results_FDR$cor > 0)] = "Higher"
table(results_FDR$direction)

## Venndiagram for intersection

load("./Analysis/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/Spearman_Relative_Genus_full_T.Rdata")
results_FDR = results[which(results$FDR < 0.1),]
results_FDR$direction = "Lower"
results_FDR$direction[which(results_FDR$cor > 0)] = "Higher"
table(results_FDR$direction)

load("./Analysis/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/Spearman_Relative_Genus_full_N.Rdata")
Nresults = results
Nresults_FDR = results[which(results$FDR < 0.1),]
Nresults_FDR$direction = "Lower"
Nresults_FDR$direction[which(Nresults_FDR$cor > 0)] = "Higher"
table(Nresults_FDR$direction)

up_T = results_FDR$Name[which(results_FDR$direction == "Higher")]
down_T = results_FDR$Name[which(results_FDR$direction == "Lower")]

up_N = Nresults_FDR$Name[which(Nresults_FDR$direction == "Higher")]
down_N = Nresults_FDR$Name[which(Nresults_FDR$direction == "Lower")]


Best_CV_coef = read.csv("./Processed_Data/Microbiome/External/9_August/Best_CV_Coefficients.csv", stringsAsFactors = FALSE)

negative_coef = Best_CV_coef$coefficients[which(Best_CV_coef$X == "blue")]
positive_coef = Best_CV_coef$coefficients[which(Best_CV_coef$X == "red")]

VennDiagram::venn.diagram(x = list(positive_coef, up_T, down_T),
                          category.names = c("positive_coef", "Tumor-distal", "Tumor-proximal"),
                          filename = paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/VennDiagram_Anatomical_location_T_",
                          "and_positive_coef.png"),
                          resolution = 600, units = "in", width = 4, height = 4)

VennDiagram::venn.diagram(x = list(negative_coef, up_T, down_T),
                          category.names = c("negative_coef", "Tumor-distal", "Tumor-proximal"),
                          filename = paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/VennDiagram_Anatomical_location_T_",
                                            "and_negative_coef.png"),
                          resolution = 600, units = "in", width = 4, height = 4)

VennDiagram::venn.diagram(x = list(positive_coef, up_N, down_N),
                          category.names = c("positive_coef", "Normal-distal", "Normal-proximal"),
                          filename = paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/VennDiagram_Anatomical_location_N_",
                                            "and_positive_coef.png"),
                          resolution = 600, units = "in", width = 4, height = 4)

VennDiagram::venn.diagram(x = list(negative_coef, up_N, down_N),
                          category.names = c("negative_coef", "Normal-distal", "Normal-proximal"),
                          filename = paste0("./Figures/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/VennDiagram_Anatomical_location_N_",
                                            "and_negative_coef.png"),
                          resolution = 600, units = "in", width = 4, height = 4)




VennDiagram::venn.diagram(x = list(up_T, down_T, up_N, down_N),
                          category.names = c("Tumor-distal", "Tumor-proximal", "Normal-distal", "Normal-proximal"),
                          filename = "./Figures/Microbiome/Re_buttal_014.2_Spearman/Relative/Genus_full/VennDiagram_Anatomical_location_TN.png",
                          resolution = 600, units = "in", width = 6, height = 4)


results_FDR_subset = results_FDR[which(results_FDR$Name %in% Best_CV_coef$covariates),]
Nresults_subset = Nresults[which(Nresults$Name %in% results_FDR_subset$Name),]

results_FDR_subset$rho_N = Nresults_subset$cor[match(results_FDR_subset$Name,
                                                     Nresults_subset$Name)]
results_FDR_subset$p_val_N = Nresults_subset$p_val[match(results_FDR_subset$Name,
                                                     Nresults_subset$Name)]
results_FDR_subset$FDR_N = Nresults_subset$FDR[match(results_FDR_subset$Name,
                                                         Nresults_subset$Name)]

results_FDR_subset$Name = gsub(".*\\D_5__", "", results_FDR_subset$Name)

results_FDR_subset= results_FDR_subset[order(results_FDR_subset$cor),]

write.csv(results_FDR_subset, file = "./Analysis/Microbiome/Re_buttal_014.2_Spearman/Overlap_with_MBR_classifier_genera_and_FDR_0.1.csv",
          row.names = FALSE)
