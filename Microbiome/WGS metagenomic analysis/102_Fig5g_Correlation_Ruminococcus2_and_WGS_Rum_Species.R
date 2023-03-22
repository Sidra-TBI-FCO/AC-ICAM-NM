
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "ComplexHeatmap", "circlize"))

# Load data
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/Relative_abundancies/clean_paired_rank_level_abundancies.Rdata")
Genus_16S = Genus_full_abundance

# Analysis
Ruminococcus_full_taxa = grep("g__Ruminococcus", rownames(Full_Taxonomy_WGS), value = TRUE)
gsub(".*s__", "", Ruminococcus_full_taxa)

Ruminococcus_names = c("s__Ruminococcus_bromii", "s__Ruminococcus_callidus", 
                       "s__Ruminococcus_champanellensis", "s__Ruminococcus_lactaris",
                       "s__Ruminococcus_sp_5_1_39BFAA")
rownames(Species_WGS)[grep("Rum", rownames(Species_WGS))]

Rum_species_WGS = Species_WGS[Ruminococcus_names,]

total_Rum_df = data.frame(Species = rownames(Rum_species_WGS), Total_reads = rowSums(Rum_species_WGS))

total_Rum_df$Percentage = total_Rum_df$Total_reads / sum(total_Rum_df$Total_reads) * 100

plot = ggplot(total_Rum_df, aes(x = "", y = Percentage, fill = Species)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "Ruminoccocus_Species", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("blue", "red", "yellow", "grey", "purple")) + coord_polar("y", start=0)

dir.create("./Figures/Microbiome/Validation_WGS_16S/102_Pie_chart_with_correlation", showWarnings = FALSE)
pdf("./Figures/Microbiome/Validation_WGS_16S/102_Pie_chart_with_correlation/102_Pie_chart_with_correlation.pdf",
    width = 10, height = 10)
plot(plot)
dev.off()

### Correlation
Rum_species_WGS_df = data.frame(t(Rum_species_WGS))
Genus_16S_df = data.frame(t(Genus_16S))

Rum_species_WGS_df$Ruminococcus2_16S = NA
Rum_species_WGS_df$Ruminococcus2_16S = Genus_16S_df$D_0__Bacteria.D_1__Firmicutes.D_2__Clostridia.D_3__Clostridiales.D_4__Ruminococcaceae.D_5__Ruminococcus.2[
  match(rownames(Rum_species_WGS_df), rownames(Genus_16S_df))
]

correlation_matrix = cor(Rum_species_WGS_df, method = "spearman")

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "darkgreen"))

pdf("./Figures/Microbiome/Validation_WGS_16S/102_Pie_chart_with_correlation/Heatmap_for_pie_colors.pdf",
    width = 6, height = 6)
Heatmap(correlation_matrix, col = col_fun)
dev.off()



