
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "ggplot2", "ggpubr", "circlize")                                                                   
ibiopak(required.bioconductor.packages)
ipak(required.packages)

# Set parameters
Tolga_pathway_excluded = "Tolga_pathway_excluded"

# define correlation significance function
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "pearson", conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
#load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Selected.pathways_ES.Rdata")
load("./Analysis/Trimmed_p/Deconvolution_and_GSEA/Selected.pathways_ES.Rdata")
ES_selected_pathways = data.frame(ES)
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
ES_Immune = t(immune_sig_df)
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
colnames(ES_selected_pathways) = gsub("X", "", colnames(ES_selected_pathways))

ES_selected_pathways = ES_selected_pathways[which(colnames(ES_selected_pathways) %in%
                                                    colnames(ES_Immune))]
colnames(ES_Immune) == colnames(ES_selected_pathways) # check if all is TRUE

ES_selected_pathways["Module11 Prolif score",] = ES_Immune["Expression Signature - Module11 Prolif score",]
ES_selected_pathways["TGFB PCA 17349583",] = ES_Immune["Expression Signature - TGFB PCA 17349583",]
ES_selected_pathways["TGFB score 21050467",] = ES_Immune["Expression Signature - TGFB score 21050467",]
ES_selected_pathways["Angiogenesis",] = ES_Immune["Expression Signature - Angiogenesis",]

exclude = c("[HM] TGF beta signaling", "[LM] Proliferation", "[HM] Angiogenesis")
ES_selected_pathways = ES_selected_pathways[-which(rownames(ES_selected_pathways) %in% exclude),]

save(ES_selected_pathways, file = "./Analysis/Trimmed_p/Deconvolution_and_GSEA/Final_Selected.pathways.Sidra.LUMC_ES.Rdata")

Rfcms$RF.predictedCMS = factor(Rfcms$RF.predictedCMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
CMSs = c("All", levels(Rfcms$RF.predictedCMS))

ES = t(data.frame(ES_selected_pathways))
rownames(ES) = gsub("\\X", "", rownames(ES))
ES = data.frame(ES)

table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in%
                                                            rownames(ES)),]
ES$ICRscore = table_cluster_assignment$ICRscore[match(rownames(table_cluster_assignment),
                                                      rownames(ES))]

ES = ES[which(rownames(ES) %in% colnames(RNASeq.QN.LOG2)),]
Rfcms = Rfcms[which(rownames(Rfcms) %in% colnames(RNASeq.QN.LOG2)),]
results = data.frame(Pathways = colnames(ES), All = NA, CMS1 = NA, CMS2 = NA, CMS3 = NA, CMS4 = NA)

cor_sign = cor.mtest(ES, 0.95)
p_mat = cor_sign[[1]]
colnames(p_mat) = colnames(ES)
rownames(p_mat) = colnames(ES)

i=1
for(i in 1:length(CMSs)){
  CMS = CMSs[i]
  if(CMS == "All"){included = rownames(Rfcms)}else{included = rownames(Rfcms)[which(Rfcms$RF.predictedCMS == CMS)]}
  ES_subset = ES[included,]
  cor_matrix = cor(ES_subset, method = "pearson")
  cor_df = data.frame(cor_matrix)
  results[, CMS] = cor_df$ICRscore[match(results$Pathways, rownames(cor_df))]

}

results$p_overall = p_mat[, "ICRscore"]
results = results[-which(results$Pathways == "ICRscore"),]

if(Tolga_pathway_excluded == "Tolga_pathway_excluded"){
  results = results[-grep("X.TPW..", results$Pathways),]
}

results = results[-grep("Breast", results$Pathways),]

results$Pathways = gsub("X.HM..", "", results$Pathways)
results$Pathways = gsub("X.IPA..", "", results$Pathways)
results$Pathways = gsub("X.TBI..", "", results$Pathways)
results$Pathways = gsub("X.TPW..", "", results$Pathways)
results$Pathways = gsub("X.LM..", "", results$Pathways)
results$Pathways = gsub("\\.", " ", results$Pathways)
results$Pathways = gsub("signaling", "sign.", results$Pathways)
results$Pathways = gsub("Signaling", "sign.", results$Pathways)
results$Pathways = gsub("Reactive oxigen species", "ROS", results$Pathways)
results$Pathways = gsub("Epithelial mesenchymal transition", "EMT", results$Pathways)
results$Pathways = gsub(" by Telomerase", "", results$Pathways)
results$Pathways = gsub(" in Eukaryotes", "", results$Pathways)

results = results[order(results$All),]

rownames(results) = results$Pathways
p_val_df = data.frame(Pathway = results$Pathways, p_val = results$p_overall, FDR = NA)
p_val_df$FDR = p.adjust(p_val_df$p_val, method = "BH", n = nrow(p_val_df))
p_val_df$p_val = -log10(p_val_df$p_val)
p_val_df$FDR = -log10(p_val_df$FDR)


results$Pathways = NULL
results$p_overall = NULL
mat = as.matrix(results)

col_fun = colorRamp2(c(min(mat), 0, max(mat)), c("#1D06FA", "white", "#00C943"))
col_fun2 = colorRamp2(c(min(p_val_df$FDR), 1.3, 2, max(p_val_df$p_val)), c("white", "yellow", "orange", "red"))
ha = HeatmapAnnotation(CMS = colnames(mat),
                       col = list(`CMS` = c("All" = "black","CMS1" = "#FF9F21", 
                                            "CMS2" = "#0074AF", "CMS3" = "#E97AA8", 
                                            "CMS4" = "#009E74")),
                       show_legend = FALSE,
                       annotation_name_gp = gpar(fontsize = 0))

ha_row = rowAnnotation(p_val = p_val_df$p_val,
                       FDR = p_val_df$FDR,
                       col = list(p_val = col_fun2,
                                  FDR = col_fun2),
                       show_legend = TRUE,
                       annotation_name_gp = gpar(fontsize = 0))

dir.create("./Figures/Trimmed_p/025_Correlation_matrix_ES", showWarnings = FALSE)
pdf(paste0("./Figures/Trimmed_p/025_Correlation_matrix_ES/025_EDF2f_v5_", 
           Tolga_pathway_excluded, "_Heatmap_Pearson_Correlation.pdf"),
    width = 5, height = 8)
Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
        col = col_fun, top_annotation = ha, show_heatmap_legend = TRUE,
        left_annotation = ha_row,
        column_title = "", row_names_max_width = unit(7, "in"),
        row_names_gp = gpar(fontsize = 12),
        show_column_names = FALSE)
dev.off()

