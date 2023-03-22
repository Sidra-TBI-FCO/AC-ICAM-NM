

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(required.packages = c("circlize", "migest", "plyr", "dplyr"))

## Set parameters
download.method = "Biolinks"

## Load data
if(download.method == "TCGA_Assembler"){
  load("./Analysis/016_CMS_Classification/TCGA_Assembler_Rfcms.Rdata")
  load("./Analysis/ICR Consensus Clustering/COAD_ICR_cluster_assignment_k2-6.Rdata")
}
if(download.method == "Biolinks"){
  load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
  load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
  load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
  table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% colnames(filtered.norm.RNAseqData)),]
}

df = data.frame(Sample = rownames(table_cluster_assignment),ICR = table_cluster_assignment$ICR_HML, CMS = NA)
df$CMS = Rfcms$RF.predictedCMS[match(df$Sample, rownames(Rfcms))]

df = df[-which(is.na(df$CMS)),]

plot_df = ddply(df, .(df$CMS, df$ICR), nrow)

colors = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
           "CMS3" = "#E97AA8", "CMS4" = "#009E74",
           "ICR High" = "red", "ICR Medium" = "green",
           "ICR Low" = "blue")

dir.create("./Figures/019_Circos_CMS_ICR", showWarnings = FALSE)
png(paste0("./Figures/019_Circos_CMS_ICR/", download.method, "_019_Circos_CMS_ICR.png"), res = 600, height = 5, width = 5, units = "in")
circos.par(gap.after = c(rep(1, length(unique(df$CMS))-1), 15, rep(1, length(unique(df$ICR))-1), 15))
chordDiagram(plot_df, order = c("CMS1", "CMS2", "CMS3", "CMS4", "ICR High", "ICR Medium", "ICR Low"), grid.col = colors,
             annotationTrack = "grid")
#circos.track(track.index = 1, panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
#}, bg.border = NA)
dev.off()
circos.clear()
