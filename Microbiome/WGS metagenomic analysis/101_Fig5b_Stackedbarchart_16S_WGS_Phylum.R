
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "dplyr", "reshape2", "colorspace", "Polychrome", "openxlsx"))

# Load data
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_paired_rank_level_abundancies.Rdata")
Phylum_16S = Phylum_abundance
rm(list=setdiff(ls(), "Phylum_16S"))
Phylum_WGS = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/WGS_validation/Microbiome_WGS_updated_167.xlsx",
                       sheet = 5)

load("./Analysis/Microbiome/Phylum_colors.Rdata")
names(Phylum_colors) = gsub("D_1__", "", names(Phylum_colors))

# Prepare data
colnames(Phylum_WGS) = gsub("\\X", "", colnames(Phylum_WGS))
rownames(Phylum_WGS) = Phylum_WGS$Phylum
Phylum_WGS$Phylum = NULL

# include only 110 in both
included = colnames(Phylum_WGS)
Phylum_16S = Phylum_16S[,included]

rownames(Phylum_16S) = gsub("D_1__", "", rownames(Phylum_16S))
rownames(Phylum_WGS) = gsub("p__", "", rownames(Phylum_WGS))

# present in both
overlap = rownames(Phylum_WGS)[which(rownames(Phylum_WGS) %in% rownames(Phylum_16S))]
unique_to_WGS = rownames(Phylum_WGS)[-which(rownames(Phylum_WGS) %in% rownames(Phylum_16S))]
unique_to_WGS

unique_to_16S = rownames(Phylum_16S)[-which(rownames(Phylum_16S) %in% overlap)]
unique_to_16S

# Only keep overlap, sum the rest and call it other
# WGS
Phylum_WGS["Other",] = NA
Phylum_WGS["Other",] = colSums(Phylum_WGS[unique_to_WGS,])
Phylum_WGS = Phylum_WGS[-which(rownames(Phylum_WGS) %in% unique_to_WGS),]

Phylum_16S["Other",] = NA
Phylum_16S["Other",] = colSums(Phylum_16S[unique_to_16S,])
Phylum_16S = Phylum_16S[-which(rownames(Phylum_16S) %in% unique_to_16S),]

# checks
colSums(Phylum_WGS)

colSums(Phylum_16S)
Phylum_16S = Phylum_16S * 100

order = colnames(Phylum_16S)[order(Phylum_16S["Fusobacteria",], decreasing = TRUE)]

Phylum_16S = Phylum_16S[, order]
Phylum_WGS = Phylum_WGS[, order]

# Generate plot
# Stacked barchart
Phylum_16S$Phylum = rownames(Phylum_16S) 
DF = melt(Phylum_16S, id = "Phylum")

plot = ggplot(DF, aes(x = variable, y = value, fill = Phylum)) + geom_bar(stat="identity") +
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

dir.create("./Figures/Exploration_reviewers/Microbiome/WGS_validation", showWarnings = FALSE)

png(paste0("./Figures/Exploration_reviewers/Microbiome/WGS_validation/Script_101_Stacked_barplot_by_Fuso_Phylum_167_samples_16S.png"),
    width = 16, height = 4, units = "in", res = 600)
plot(plot)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/WGS_validation/Script_101_Stacked_barplot_by_Fuso_Phylum_167_samples_16S.pdf"),
    width = 16, height = 4)
plot(plot)
dev.off()


## WGS
Phylum_WGS$Phylum = rownames(Phylum_WGS)
DF2 = melt(Phylum_WGS, id = "Phylum")

plot = ggplot(DF2, aes(x = variable, y = value, fill = Phylum)) + geom_bar(stat="identity") +
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

dir.create("./Figures/Exploration_reviewers/Microbiome/WGS_validation", showWarnings = FALSE)

png(paste0("./Figures/Exploration_reviewers/Microbiome/WGS_validation/Script_101_Stacked_barplot_by_Fuso_Phylum_167_samples_WGS.png"),
    width = 16, height = 4, units = "in", res = 600)
plot(plot)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/WGS_validation/Script_101_Stacked_barplot_by_Fuso_Phylum_167_samples_WGS.pdf"),
    width = 16, height = 4)
plot(plot)
dev.off()




