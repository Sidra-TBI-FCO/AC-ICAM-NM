

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggpubr"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")

length(unique(finalMafFiltered$Patient_ID))
table(finalMafFiltered$Variant_Classification)

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
# rownames(Genus_full_abundance)[grep("Fuso", rownames(Genus_full_abundance))
Genus = "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
Top10 = "Top10" # "Top10"


dir.create("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis", showWarnings = FALSE)
dir.create(paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status"), showWarnings = FALSE)
dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_06_by_specific_mutation", showWarnings = FALSE)

patients = substring(colnames(Genus_full_abundance), 1, 3)
finalMafFiltered = finalMafFiltered[which(finalMafFiltered$Patient_ID %in% patients),]

length(unique(finalMafFiltered$Patient_ID))
if(Top10 == ""){
  tbl = data.frame(table(finalMafFiltered$Hugo_Symbol))
  genes = as.character(tbl$Var1[which(tbl$Freq > 10)])
}else{
  load(paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status/results_df_",
              Genus, ".Rdata"))
  results_df = results_df[order(results_df$p_value_Wilcoxon, decreasing = FALSE),]
  genes = results_df$Gene[1:10]
}


df = data.frame(Patient_ID = substring(colnames(Genus_full_abundance), 1, 3),
                Abundance = Genus_full_abundance[Genus,])

results_df = data.frame(Gene = genes, p_value_Wilcoxon = NA)

i=1
#for (i in 1:length(genes)){
  #gene = genes[i]
  gene = "KMT2C"
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
  df_plot$Microbiome_Category = "Positive"
  df_plot$Microbiome_Category[which(df_plot$Abundance == 0)] = "Negative"
  
  table(df_plot$Mut, df_plot$Microbiome_Category)
  chisq.test(table(df_plot$Mut, df_plot$Microbiome_Category))
  
  # Prepare data
  DF1 <- df_plot %>%
    group_by(Mut, Microbiome_Category) %>%
    summarise(count=n()) %>%
    mutate(perc=count/sum(count))
  
  
  plot = ggplot(DF1, aes(x =  Mut, y =perc*100, fill = Microbiome_Category)) + geom_bar(stat="identity") +
    labs(x = "", y = "Percentage", fill = "", face = "bold") +
    theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
          panel.background = element_rect(fill = "white", colour = "grey"),
          axis.text = element_text(size = 19, colour = "black"),
          axis.title = element_text(size = 19, colour = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 19, colour = "black"),
          legend.position = "none")+
    scale_fill_manual(values= c("Positive" = "#FF4500", "Negative" = "#6495ED"))
    
    png(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_06_by_specific_mutation/Stacked_barchart_", Genus, 
               "_by_", gene, "_mutation.png"), res = 600, width = 3, height = 4, units = "in")
    plot(plot)
    dev.off()
  }
  


save(results_df, file = paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/MOI_06_by_specific_mutation_status/results_df_", Genus,".Rdata"))



