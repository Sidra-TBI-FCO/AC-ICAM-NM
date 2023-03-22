
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(required.packages = c("ggplot2", "dplyr"))

# Set parameters
exclude_medium = "exclude_medium" # "include_medium" or "exclude_medium" or "include_medium_cat"
stage_filter = ""

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load(paste0("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/frequency_df.Rdata"))
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

Merged_dataset$ICRscore = table_cluster_assignment$ICRscore[match(Merged_dataset$Patient_ID,
                                                                  substring(rownames(table_cluster_assignment), 1, 3))]
Merged_dataset$ICR_cluster = table_cluster_assignment$ICR_HML[match(Merged_dataset$Patient_ID,
                                                                    substring(rownames(table_cluster_assignment), 1, 3))]

Merged_dataset$CMS = Rfcms$RF.predictedCMS[match(Merged_dataset$Patient_ID, substring(rownames(Rfcms), 1, 3))]

dim(RNASeq.QN.LOG2)
table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(RNASeq.QN.LOG2)),]
median_ICRscore = median(table_cluster_assignment$ICRscore)

if(exclude_medium == "exclude_medium"){
  Merged_dataset = Merged_dataset[-which(Merged_dataset$ICR_cluster == "ICR Medium"),]
  color_vals = c("ICR High" = "red", "ICR Low" = "blue")
}

if(exclude_medium == "include_medium"){
  Merged_dataset$ICR_by_median = NA
  Merged_dataset$ICR_by_median[which(Merged_dataset$ICRscore < median_ICRscore)] = "ICR Low"
  Merged_dataset$ICR_by_median[which(Merged_dataset$ICRscore >= median_ICRscore)] = "ICR High"
  table(Merged_dataset$ICR_by_median, Merged_dataset$Immunoedited)
  Merged_dataset$ICR_cluster = Merged_dataset$ICR_by_median
  color_vals = c("ICR High" = "red", "ICR Low" = "blue")
}

table(Merged_dataset$ajcc_pathologic_tumor_stage, exclude = NULL)
if(stage_filter == ""){}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage == stage_filter),]
}

TCR_Overview = TCR_Overview[which(TCR_Overview$Tissue == "Primary tumor"),]
Merged_dataset$TCR_productive_clonality = TCR_Overview$productive_clonality[match(Merged_dataset$Patient_ID,
                                                                                  TCR_Overview$Patient_ID)]

Merged_dataset$Immunoedited_ICR = paste(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)

# Scatterplot
scatterplot = ggplot(data = Merged_dataset, aes(x = Immunoediting_score, y = ICRscore)) +
  geom_point(aes(color = Immunoedited_ICR)) +
  scale_color_manual(values = c("immunoedited ICR High" = "#FF3806",
                                "less immunoedited ICR High" = "#FF9F00",
                                "immunoedited ICR Low" = "#009CFC",
                                "less immunoedited ICR Low" = "#5233FC")) +
  #stat_cor(method = "pearson", size = 6) +
  #scale_y_log10() +
  theme_bw() +
  ylab("ICR score") +
  xlab("Immunoediting score") +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        aspect.ratio = 1/1) #+

pdf(paste0("./Figures/WES/022_Ratio_Neoantigens_Mutational_load_exploration/Fig4c_Scatterplot_ICRscore_by_immunoediting_status_",
           exclude_medium, stage_filter, ".pdf"),
    height = 5, width = 5)
plot(scatterplot)
dev.off()

