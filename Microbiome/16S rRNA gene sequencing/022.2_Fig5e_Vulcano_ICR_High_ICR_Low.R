
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggrepel"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
subset = ""
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)

# Load data
load(paste0("./Analysis/Microbiome/008.2_ICR_wilcoxon_boxplot/Relative/Genus_full/New_filter_138_genera__by_ICR_cluster_Relative_Genus_full.Rdata"))

# Make Volcano plot
results = results[-which(results$Name == "Unknown"),]

results$FDR = p.adjust(results$p_val, method = "fdr", n = nrow(results))
results = results[order(results$p_val),]

results$Direction = NA
results$Direction[which(results$mean_ICR_Low > results$mean_ICR_High)] = "Enriched in ICR Low"
results$Direction[which(results$mean_ICR_Low <= results$mean_ICR_High)] = "Enriched in ICR High"
#results$Ratio_mean = results$mean_ICR_High / results$mean_ICR_Low

results$Name = as.character(results$Name)
results$p_value_log10 = -log10(results$p_val)
results$FDR_log10 = -log10(results$FDR)
results$color = "non-significant"
results$color[results$FDR < 0.1 & results$Direction == "Enriched in ICR Low"] = "p < 0.05 enriched in ICR Low"
results$color[results$FDR < 0.1 & results$Direction == "Enriched in ICR High"] = "p < 0.05 enriched in ICR High"

#results = results[order(results$ratio, decreasing = TRUE),]
results$label = gsub(".*D_5__", "", results$Name)
results$Phylum = gsub(".*D_1__", "", results$Name)
results$Phylum = gsub("\\ D_2__.*", "", results$Phylum)
#results$label = paste(results$label, " (", results$Phylum, ")", sep = "")
results$Plot_label = results$label
results$Plot_label[-which(results$FDR < 0.1)] = NA
#results$Plot_label = NA

#results$Plot_label[which(results$p_value_log10 < 10)] = NA
#results$Plot_label[which(results$log_ratio > 0.01)] = results$label[which(results$log_ratio > 0.01)]


#results = results[-which(results$mean_N < 0.000001 & results$mean_T < 0.00001),] # delete low abundance rank level microbiota

results$Ratio_mean_log10 = log(results$Ratio_mean,10)
#results = results[-which(results$Ratio_mean_log10 == Inf),]
#results = results[-which(results$Ratio_mean_log10 == -Inf),]
#results$Plot_label[which(results$Ratio_mean_log10 < 0.5 & results$p_value_log10 < 10)] = NA

plotA = ggplot(results, aes(x= FC_log2, y= FDR_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "p < 0.05 enriched in ICR Low" = "blue", "p < 0.05 enriched in ICR High" = "red")) +
  theme_bw() +
  #xlim(-2, 2) +
  ylab("-log10 FDR") +
  xlab("") +
  geom_line(y = 1, linetype = "dashed", color = "darkgrey") +
  #geom_line(y = 2.1, linetype = "dashed", color = "darkgrey") + # change this for final version (FDR = 0.05)
  geom_text_repel(aes(x = FC_log2, FDR_log10, 
                      label = Plot_label), size = 3,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none",
        aspect.ratio = 1/1.5)

dir.create("./Figures/Microbiome/008.3_ICR_High_Low_Volcano", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/008.3_ICR_High_Low_Volcano/FDR_New_Filter_138_Genera_No_label_ICR_High_Low_Volcano_mean_ratio_", Type, "_abundance_", Rank, ".png"), 
    width = 5, height = 4, units = "in", res = 600)
plot(plotA)
dev.off()

pdf(paste0("./Figures/Microbiome/008.3_ICR_High_Low_Volcano/FDR_New_Filter_138_Genera_No_label_ICR_High_Low_Volcano_mean_ratio_", Type, "_abundance_", Rank, ".pdf"), 
    width = 5, height = 4)
plot(plotA)
dev.off()




