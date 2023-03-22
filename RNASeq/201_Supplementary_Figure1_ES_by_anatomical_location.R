
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

dir.create("./Figures/Exploration_reviewers/ES_by_anatomical_location", showWarnings = FALSE)
dir.create("./Analysis/Exploration_reviewers/ES_by_anatomical_location", showWarnings = FALSE)

ipak(c("RColorBrewer", "ggplot2", "easyGgplot2"))

# Set parameters
Gene_set = "ConsensusTME"
Include_CMS = ""

# Load data
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v2_All_Immune_gene_signatures_table.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")

if(Gene_set == "ConsensusTME"){
  Tfh = immune_sig_df[, "Bindea |  TFH", drop = FALSE]
  load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/ConsensusTME_COAD_ES.Rdata")
  immune_sig_df = t(ES) 
  immune_sig_df = cbind(immune_sig_df, Tfh)
  immune_sig_df$ICRscore = table_cluster_assignment$ICRscore[match(rownames(immune_sig_df), rownames(table_cluster_assignment))]
  immune_sig_df = immune_sig_df[which(substring(rownames(immune_sig_df), 1, 3) %in% clinical_data$Patient_ID),]
}

if(Include_CMS == ""){}else{
  clinical_data$CMS = Rfcms$RF.predictedCMS[match(clinical_data$Patient_ID, substring(rownames(Rfcms), 1, 3))]
  clinical_data = clinical_data[which(clinical_data$CMS == Include_CMS),]
}

# Spearman correlation anatomical location and enrichment score
df = data.frame(Patient_ID = clinical_data$Patient_ID, Anatomic_location = clinical_data$tumour_anatomic_site,
                ES = NA)
df$Anatomic_location[which(df$Anatomic_location == "ceceum")] = "caecum"

df$Anatomic_location = factor(df$Anatomic_location, levels = c("caecum", "colon ascendens", "flexura hepatica", 
                                                                   "colon transversum", "flexura lienalis", "colon descendens", 
                                                                   "colon sigmoideum", "rectosigmoideum"))
table(df$Anatomic_location)
df$Anatomic_location = as.numeric(df$Anatomic_location)

results_df = data.frame(Pathway = colnames(immune_sig_df), Spearman_cor = NA, p_value = NA)

i= 1 #45
for (i in 1:ncol(immune_sig_df)){
  sign = colnames(immune_sig_df)[i]
  df$ES = immune_sig_df[, sign][match(df$Patient_ID, substring(rownames(immune_sig_df), 1, 3))]
  cor = cor.test(df$ES, df$Anatomic_location, method = "spearman")
  results_df$Spearman_cor[which(results_df$Pathway == sign)] = cor$estimate
  results_df$p_value[which(results_df$Pathway == sign)] = cor$p.value
  
  df_plot = df
  df_plot$Anatomic_location = factor(df_plot$Anatomic_location)
  levels(df_plot$Anatomic_location) = c("caecum", "colon ascendens", "flexura hepatica", 
                      "colon transversum", "flexura lienalis", "colon descendens", 
                      "colon sigmoideum", "rectosigmoideum")
  
  plot = ggplot(df_plot, aes(x = Anatomic_location, y = ES, fill = Anatomic_location)) +
    scale_fill_manual(values = brewer.pal(n = 8, name = "Set3")) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 0.8) +
    ylab(gsub("_", " ", sign)) +
    theme_bw() +
    xlab("") +
    theme(axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 15, colour = "black",
                                     angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),
          legend.position = "none") +
    geom_text(x = 5, y = max(df_plot$ES), label = paste0("Rho = ", round(cor$estimate, 2), ", ", "p  = ", 
                                                         formatC(cor$p.value, format = "e", digits = 2)),
              check_overlap = TRUE)
  
  sign = gsub("\\/", "-", sign)
  
  assign(paste0("p", i, sep = ""), plot)
  png(paste0("./Figures/Exploration_reviewers/ES_by_anatomical_location/", Include_CMS, "_Boxplot_ES_", sign, "_by_anatomic_location.png"),
      res = 600, units = "in", width = 5, height = 4)
  plot(plot)
  dev.off()
  
}

plots = paste("p", 1:ncol(immune_sig_df), sep = "")
list_of_plots = mget(plots)

results_df$FDR = p.adjust(results_df$p_value, method = "BH", n = nrow(results_df))

list_of_plots = list_of_plots[order(results_df$Spearman_cor)][1:20]
results_df = results_df[order(results_df$Spearman_cor),]


pdf(paste0("./Figures/Exploration_reviewers/ES_by_anatomical_location/New_", Include_CMS,"_Boxplot_ES_", sign, "_",
           Gene_set,
           ".pdf"),
    width = 17, height = 23)
ggplot2.multiplot(plotlist = list_of_plots[1:ncol(immune_sig_df)], cols=4)
dev.off()

save(results_df, file = paste0("./Analysis/Exploration_reviewers/ES_by_anatomical_location/Spearman_cor_results_df_", Gene_set,".Rdata"))
write.csv(results_df, file = paste0("./Analysis/Exploration_reviewers/ES_by_anatomical_location/", Include_CMS, "_Spearman_cor_results_df_", Gene_set, ".csv"),
          row.names = FALSE)

