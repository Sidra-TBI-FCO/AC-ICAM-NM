

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ggrepel"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full"
# for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_Wilcoxon_test/48_TN_Paired_Wilcoxon_test.Rdata")

# Make Volcano plot
results$Name = as.character(results$Name)
results$p_value_log10 = -log10(results$p_val)
results$FDR_log10 = -log10(results$FDR)
results$color = "non-significant"
results$color[results$FDR < 0.05 & results$Direction == "Enriched in tumor"] = "FDR < 0.1 enriched in Tumor"
results$color[results$FDR < 0.05 & results$Direction == "Enriched in normal"] = "FDR < 0.1 enriched in Normal"
results$Ratio_mean_log10 = log(results$Ratio_mean,10)
results$Ratio_mean_log2 = log(results$Ratio_mean,2)

#results = results[order(results$ratio, decreasing = TRUE),]
results$label = gsub(".*D_5__", "", results$Name)
results$Phylum = gsub(".*D_1__", "", results$Name)
results$Phylum = gsub("\\ D_2__.*", "", results$Phylum)
#results$label = paste(results$label, " (", results$Phylum, ")", sep = "")
results$Plot_label = results$label
results$Plot_label[-which(results$FDR < 0.05)] = NA
#results$Plot_label[which(results$Ratio_mean_log2 <= 1)] = NA

#results$Plot_label[which(results$log_ratio_paired < 0.5)] = NA
#results$Plot_label[which(results$p_value_log10 < 10)] = NA
#results$Plot_label[which(results$log_ratio > 0.01)] = results$label[which(results$log_ratio > 0.01)]

results = results[-which(results$Name == "Unknown"),]

#results$Plot_label[which(results$Ratio_mean_log2 < 0.5 & results$FDR_log10 < 1.30103)] = NA

plotB = ggplot(results, aes(x= Ratio_mean_log10, y= p_value_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  xlim(-2, 2) +
  ylab("-log10 pvalue") +
  xlab("Log10 ratio \n (mean tumor / mean normal)") +
  geom_line(y = 1.3, linetype = "dashed", color = "darkgrey") +
  #geom_line(y = 2.1, linetype = "dashed", color = "darkgrey") + # change this for final version (FDR = 0.05)
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = 0.5, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = -0.5, linetype = "dashed", color = "darkgrey") +
  geom_text_repel(aes(x = Ratio_mean_log10, p_value_log10, 
                      label = Plot_label), size = 0,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/48_Paired_wilcoxon_test_TN_Volcano_plot.png"), 
    width = 5, height = 4, units = "in", res = 600)
plot(plotB)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/48_Paired_wilcoxon_test_TN_Volcano_plot.pdf"), width = 5, height = 4)
plot(plotB)
dev.off()

plotC = ggplot(results, aes(x= Ratio_mean_log2, y= FDR_log10, colour= color, label = Plot_label)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("non-significant" = "black", "FDR < 0.1 enriched in Tumor" = "#1E90FF", "FDR < 0.1 enriched in Normal" = "#FF00FF")) +
  theme_bw() +
  ylim(0, 20) +
  #xlim(-2, 2) +
  ylab("-log10 FDR") +
  xlab("Log2 ratio \n (mean tumor / mean normal)") +
  geom_line(y = 1.30103, linetype = "dashed", color = "darkgrey") +
  #geom_line(y = 2.1, linetype = "dashed", color = "darkgrey") + # change this for final version (FDR = 0.05)
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = 0.5, linetype = "dashed", color = "darkgrey") +
  #geom_line(x = -0.5, linetype = "dashed", color = "darkgrey") +
  geom_text_repel(aes(x = Ratio_mean_log2, y = FDR_log10, 
                      label = Plot_label), size = 3, max.overlaps = 15,
                  segment.size = 0.2) +
  theme(axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        legend.position = "none")

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/FDR_48_Paired_wilcoxon_test_TN_Volcano_plot.png"), 
    width = 5, height = 4, units = "in", res = 600)
plot(plotC)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/48_Paired_wilcoxon_test_TN/FDR_48_Paired_wilcoxon_test_TN_Volcano_plot.pdf"), width = 5, height = 4)
plot(plotC)
dev.off()
