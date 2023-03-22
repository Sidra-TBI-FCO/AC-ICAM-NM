# Vulcano plot

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "ggrepel")
ipak(required.packages)

set.seed(4)

# Set parameters
signif = "FDR" # "FDR" or "p-value"

# load data
results = read.csv("./Analysis/TCR/009_Correlation_gene_expression_clonality/114_pearson_correlation_genes_productive_clonality.csv",
                   stringsAsFactors = FALSE)
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))

# Analysis
dir.create("./Figures/TCR/010_Vulcano_plots", showWarnings = FALSE)

results$FDR = p.adjust(results$p_value, method = "BH", n = nrow(results))
if(signif == "FDR"){
  results$p_value = results$FDR
}

results$p_value_log10 = -log10(results$p_value)
results$color = "non-significant"
results$color[results$p_value < 0.05 & results$Corr_coefficient > 0] = "p < 0.05_pos"
results$color[results$p_value < 0.05 & results$Corr_coefficient < 0] = "p < 0.05_neg"


results = results[order(results$Corr_coefficient, decreasing = TRUE),]
results$Gene_label = NA
results$Gene_label[1:10] = results$Gene[1:10]
top20 = results$Gene_label[1:10]
intersect(top20, ICR_genes)

results = results[order(results$Corr_coefficient),]
results$Gene_label[1:10] = results$Gene[1:10]

results$Gene[which(results$Gene == "TNFRSF9")] = "CD137"
results$Gene[which(results$Gene == "ENTPD1")] = "CD39"
results$Gene[which(results$Gene == "ITGAE")] = "CD103"
results$Gene[which(results$Gene == "TNFRSF18")] = "GITR"
results$Gene[which(results$Gene == "HAVCR2")] = "TIM3"

selected_genes = c("LAG3", "PDCD1", "TIM3",
                   "GITR", "CD103", "CD137", "CD39", "CXCL13")
#results$Gene_label[which(results$Gene %in% selected_genes)] = results$Gene[which(results$Gene %in% selected_genes)]

results_sub = results[which(results$Gene %in% c(ICR_genes, selected_genes)), ]
dir.create("./Analysis/TCR/010_TCR_ICR_tumor_reactive_cells_subset", showWarnings = FALSE)
write.csv(results_sub, file = "./Analysis/TCR/010_TCR_ICR_tumor_reactive_cells_subset/results_sub.csv",
          row.names = FALSE)


results$color2 = results$color
results$color2[which(results$Gene %in% selected_genes)] = "selected_genes"
results$color2[which(results$Gene %in% ICR_genes)] = "ICR"

plot = ggplot(results, aes(x= Corr_coefficient, y= p_value_log10, colour= color, label = Gene_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "ICR" = "orange", "p < 0.05_neg" = "#5287DB",
                                "selected_genes" = "#00B050", "p < 0.05_pos" = "red")) +
  theme_bw() +
  xlim(-0.6, 0.702) +
  ylab(paste0("log10 ", signif)) +
  xlab("Pearson's r") +
  geom_text_repel(aes(x = Corr_coefficient, p_value_log10, 
                      label = Gene_label, color = color), size = 3.6,
                  max.overlaps = 30) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

#png("./Figures/TCR/010_Vulcano_plots/v4_114_Vulcano_p_0.001.png", width = 5, height = 4, units = "in", res = 600)
#plot(plot)
#dev.off()

pdf(paste0("./Figures/TCR/010_Vulcano_plots/Fig2d_114_Vulcano_p_0.05_", signif, ".pdf"), width = 5, height = 4)
plot(plot)
dev.off()
