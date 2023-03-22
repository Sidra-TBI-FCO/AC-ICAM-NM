

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
Fusobacterium_full_taxa = grep("g__Fusobacterium", rownames(Full_Taxonomy_WGS), value = TRUE)
gsub(".*s__", "", Fusobacterium_full_taxa)

Fuso_names = c("s__Fusobacterium_gonidiaformans", "s__Fusobacterium_mortiferum", 
                       "s__Fusobacterium_necrophorum", "s__Fusobacterium_nucleatum",
                       "s__Fusobacterium_periodonticum", "s__Fusobacterium_ulcerans",
                       "s__Fusobacterium_varium")


rownames(Species_WGS)[grep("Fuso", rownames(Species_WGS))]

Fuso_species_WGS = Species_WGS[Fuso_names,]

total_Fuso_df = data.frame(Species = rownames(Fuso_species_WGS), Total_reads = rowSums(Fuso_species_WGS))

total_Fuso_df$Percentage = total_Fuso_df$Total_reads / sum(total_Fuso_df$Total_reads) * 100

plot = ggplot(total_Fuso_df, aes(x = "", y = Percentage, fill = Species)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "Species", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("blue", "red", "yellow", "grey", "purple", "green", "orange")) + coord_polar("y", start=0)

dir.create("./Figures/Microbiome/Validation_WGS_16S/103_Pie_chart_with_correlation", showWarnings = FALSE)
pdf("./Figures/Microbiome/Validation_WGS_16S/103_Pie_chart_with_correlation/103_Pie_chart_with_correlation_Fusobacterium.pdf",
    width = 10, height = 10)
plot(plot)
dev.off()

### Correlation
Fuso_species_WGS_df = data.frame(t(Fuso_species_WGS))
Genus_16S_df = data.frame(t(Genus_16S))

Fuso_species_WGS_df$Fuso_16S = NA
Fuso_species_WGS_df$Fuso_16S = Genus_16S_df$D_0__Bacteria.D_1__Fusobacteria.D_2__Fusobacteriia.D_3__Fusobacteriales.D_4__Fusobacteriaceae.D_5__Fusobacterium[
  match(rownames(Fuso_species_WGS_df), rownames(Genus_16S_df))
]

correlation_matrix = cor(Fuso_species_WGS_df, method = "spearman")

cor_test = cor.test(Fuso_species_WGS_df$s__Fusobacterium_nucleatum, Fuso_species_WGS_df$Fuso_16S, method = "spearman")
cor_test$p.value

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "darkorange"))

pdf("./Figures/Microbiome/Validation_WGS_16S/103_Pie_chart_with_correlation/Heatmap_for_pie_colors_Fuso.pdf",
    width = 6, height = 6)
Heatmap(correlation_matrix, col = col_fun)
dev.off()



