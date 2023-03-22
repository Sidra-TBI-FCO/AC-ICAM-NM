

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(required.packages = c("circlize", "plyr", "dplyr", "stringr")) #"migest"

## Set parameters
exclude = c("Conpair_lower_90_percent", "non-epithelial")

## Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

#Analysis
table_cluster_assignment = table_cluster_assignment[-which(substring(rownames(table_cluster_assignment), 1, 3) %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]),]
df = data.frame(Sample = rownames(table_cluster_assignment),ICR = table_cluster_assignment$ICR_HML, CMS = NA)
df$CMS = Rfcms$RF.predictedCMS[match(df$Sample, rownames(Rfcms))]

df = df[-which(is.na(df$CMS)),]

plot_df = ddply(df, .(df$CMS, df$ICR), nrow)

colors = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
           "CMS3" = "#E97AA8", "CMS4" = "#009E74",
           "ICR High" = "red", "ICR Medium" = "green",
           "ICR Low" = "blue")

dir.create("./Figures/Trimmed_p/019_Circos_CMS_ICR", showWarnings = FALSE)
pdf("./Figures/Trimmed_p/019_Circos_CMS_ICR/019_Circos_CMS_ICR.pdf", height = 5, width = 5)
circos.par(gap.after = c(rep(1, length(unique(df$CMS))-1), 15, rep(1, length(unique(df$ICR))-1), 15))
chordDiagram(plot_df, order = c("CMS1", "CMS2", "CMS3", "CMS4", "ICR High", "ICR Medium", "ICR Low"), grid.col = colors,
             annotationTrack = "grid")
#circos.track(track.index = 1, panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#              facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
#}, bg.border = NA)
dev.off()
circos.clear()
