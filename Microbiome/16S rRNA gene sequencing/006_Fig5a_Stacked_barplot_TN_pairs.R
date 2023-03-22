
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "dplyr", "reshape2", "colorspace", "Polychrome"))

# Load data
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_paired_rank_level_abundancies.Rdata")

# Set parameters
Tissue = "T"
ordering = "by_delta_Fuso"

# Stacked barchart
Phylum_abundance$Phylum = rownames(Phylum_abundance) 
DF = melt(Phylum_abundance, id = "Phylum")

Phylum_abundance = as.matrix(Phylum_abundance)
mode(Phylum_abundance) = "numeric"
abundance_overall = data.frame(rowSums(Phylum_abundance, na.rm = TRUE))
abundance_overall$color = NA
abundance_overall$color[which(abundance_overall$rowSums.Phylum_abundance..na.rm...TRUE. < 0.1)] = "darkgrey"
delete = rownames(abundance_overall)[which(abundance_overall$color == "darkgrey")]

data("alphabet")
names(alphabet)[c(24:26, 1:23)] = unique(DF$Phylum)
alphabet[which(names(alphabet) == "D_1__Verrucomicrobia")] = "#2ED9FF"
alphabet[which(names(alphabet) == "D_1__Nanoarchaeaeota")] = "#FC1CBF"  
alphabet[which(names(alphabet) == "D_1__Fusobacteria")] = "#FEAF16"
alphabet[which(names(alphabet) == "D_1__Spirochaetes")] = "#FC1CBF"
alphabet[which(names(alphabet) %in% delete)] = "darkgrey"

Phylum_colors = alphabet

save(Phylum_colors, file = "./Analysis/Microbiome/Phylum_colors.Rdata")

DF$Tissue = substring(DF$variable, 4,4)

if(ordering == "by_delta_Fuso"){
  DF_T = DF[which(DF$Tissue == "T" & DF$Phylum == "D_1__Fusobacteria"),]
  DF_N = DF[which(DF$Tissue == "N" & DF$Phylum == "D_1__Fusobacteria"),]
  new_paired = DF_T
  new_paired$variable = as.character(new_paired$variable)
  DF_N$variable = as.character(DF_N$variable)
  new_paired$value_in_N = DF_N$value[match(substring(new_paired$variable, 1, 3), substring(DF_N$variable, 1, 3))]
  new_paired$variable_N = DF_N$variable[match(substring(new_paired$variable, 1, 3), substring(DF_N$variable, 1, 3))]
  new_paired$Delta = new_paired$value - new_paired$value_in_N
  new_paired = new_paired[order(new_paired$Delta, decreasing = TRUE),]
  if(Tissue == "T"){
    order = new_paired$variable
  }
  if(Tissue == "N"){
    order = new_paired$variable_N
  }
}

if(ordering == "by_ICR"){
  load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
  DF2 = DF
  DF2$Patient = substring(DF$variable, 1, 3)
  DF2$ICR = table_cluster_assignment$ICRscore[match(DF2$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
  DF2$ICR_cluster = table_cluster_assignment$ICR_HML[match(DF2$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
  DF2 = DF2[order(DF2$ICR),]
  DF2$ICR_cluster = factor(DF2$ICR_cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))
  DF2 = DF2[order(DF2$ICR_cluster),]
  DF2 = DF2[which(DF2$Tissue == Tissue),]
  order = DF2$variable
}

#DF = DF[order(DF$Tissue),]
DF = DF[which(DF$Tissue == Tissue),]

if(ordering == "by_delta_Fuso"){
  DF$variable = as.character(DF$variable)
  DF$variable = factor(DF$variable, levels = order)
}else{
  DF$variable = factor(DF$variable, levels = unique(DF$variable))
}

DF = DF[order(DF$variable),]

plot = ggplot(DF, aes(x = variable, y = value*100, fill = Phylum)) + geom_bar(stat="identity") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_fill_manual(values = Phylum_colors) +
  ylab("Percent") +
  xlab("Sample")
  

dir.create("./Figures/Microbiome", showWarnings = FALSE)
dir.create("./Figures/Microbiome/006_Stacked_barplot_TN_pairs", showWarnings = FALSE)

png(paste0("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/v3_", ordering, Tissue, "_Phylum_all_samples.png"),
    width = 16, height = 4, units = "in", res = 600)
plot(plot)
dev.off()

pdf(paste0("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/v3_", ordering, Tissue, "_Phylum_all_samples.pdf"),
    width = 16, height = 4)
plot(plot)
dev.off()

DF = DF[-which(DF$Phylum %in% delete),]

legend_plot = ggplot(DF, aes(x = variable, y = value*100, fill = Phylum)) + geom_bar(stat="identity") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = Phylum_colors) +
  ylab("Percent") +
  xlab("Sample")

png("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/Phylum_legend.png",
    width = 8, height = 4, units = "in", res = 600)
plot(legend_plot)
dev.off()

pdf("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/Phylum_legend.pdf",
    width = 8, height = 4)
plot(legend_plot)
dev.off()

## Specific samples graph
DF = DF[which(DF$variable %in% c("002N", "002T", "009N", "009T", "010N", "010T", "020N", "020T", "022N", "022T", "023N", "023T")),]

plot = ggplot(DF, aes(x = variable, y = value*100, fill = Phylum)) + geom_bar(stat="identity") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black",
                                   angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_fill_manual(values = Phylum_colors) +
  ylab("Percent") +
  xlab("Sample")

png("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/v2_Phylum_specific_samples.png",
    width = 8, height = 4, units = "in", res = 600)
plot(plot)
dev.off()

pdf("./Figures/Microbiome/006_Stacked_barplot_TN_pairs/v2_Phylum_specific_samples.pdf",
    width = 8, height = 4)
plot(plot)
dev.off()


