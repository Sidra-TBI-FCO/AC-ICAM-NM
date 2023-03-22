

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")

length(unique(finalMafFiltered$Patient_ID))
table(finalMafFiltered$Variant_Classification)

# Set parameters 
rownames(Species_WGS)
Species = "s__Fusobacterium_nucleatum"
Top10 = "Top10" # "Top10"

dir.create("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status"), showWarnings = FALSE)
dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_06_by_specific_mutation", showWarnings = FALSE)

# Analysis
Species_WGS = as.matrix(Species_WGS)
Species_WGS = Species_WGS / 100

patients = substring(colnames(Species_WGS), 1, 3)
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% patients),]

length(unique(finalMafFiltered$Patient_ID))
if(Top10 == ""){
  tbl = data.frame(table(finalMafFiltered$Hugo_Symbol))
  genes = as.character(tbl$Var1[which(tbl$Freq > 10)])
}else{
  load(paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status/WGS_results_df_",
              Species, ".Rdata"))
  results_df = results_df[order(results_df$p_value_Wilcoxon, decreasing = FALSE),]
  genes = results_df$Gene[1:10]
}


df = data.frame(Patient_ID = substring(colnames(Species_WGS), 1, 3),
                Abundance = Species_WGS[Species,])

results_df = data.frame(Gene = genes, p_value_Wilcoxon = NA)

i=1
for (i in 1:length(genes)){
  #gene = genes[i]
  gene = "BRAF"
  subset = finalMafFiltered[which(finalMafFiltered$Hugo_Symbol == gene),]
  mut_samples = unique(subset$Patient_ID)
  df$Mut = "No mutation"
  df$Mut[which(df$Patient_ID %in% mut_samples)] = "Mutation"
  df$Mut = factor(df$Mut, levels = c("No mutation", "Mutation"))
  
  wilcox = wilcox.test(df$Abundance[which(df$Mut == "No mutation")], df$Abundance[which(df$Mut == "Mutation")])
  results_df$p_value_Wilcoxon[which(results_df$Gene == gene)] = wilcox$p.value
  results_df$Mean_mut[which(results_df$Gene == gene)] = mean(df$Abundance[which(df$Mut == "Mutation")])
  results_df$Mean_no_mut[which(results_df$Gene == gene)] = mean(df$Abundance[which(df$Mut == "No mutation")])
  
  if(length(mut_samples) < 10){next}
  df_plot = df
  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
  
  
  if(Top10 == "Top10"){
    plot = ggplot(df_plot, aes(x = Mut, y = Abundance)) +
      geom_boxplot(outlier.shape = NA, aes(fill=Mut)) +
      scale_fill_manual(values = c("No mutation" = "#51EE77", "Mutation" = "#FD76FB")) +
      geom_jitter(width = 0.2, size = 0.8) +
      theme_bw() +
      stat_compare_means() +
      scale_y_log10() +
      ylab("") +
      xlab("") +
      theme(axis.text.y = element_text(size = 15, colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 15, colour = "black"),
            axis.title.y = element_text(size = 15, colour = "black"),
            legend.position = "none")
    
    pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_06_by_specific_mutation/EDF8l_WGS_based_abundance_", Species, 
               "_by_", gene, "_mutation.pdf"), width = 3, height = 4)
    plot(plot)
    dev.off()
  }
  
}

save(results_df, file = paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status/WGS_results_df_", Species,".Rdata"))



