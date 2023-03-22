#################################################################
###
### Create RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA
###
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
Gene_set = "DDR" # "Selected.pathways"

# Create folders and log file
dir.create("./Figures",showWarnings = FALSE)
dir.create("./Figures/Heatmaps", showWarnings = FALSE)
dir.create(paste0("./Figures/Heatmaps/", Gene_set, "_Heatmaps"), showWarnings = FALSE)
dir.create("./Analysis/Deconvolution_and_GSEA", showWarnings = FALSE)

#load data
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
load(paste0(toolbox.path, "/GSEA list/immune.gene.lists.v4.Rdata"))
load(paste0(toolbox.path, "/Selected Pathways/Selected.pathways.3.4.RData"))
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))
load(paste0(toolbox.path, "/GSEA list/DDR.gene.list.Rdata"))

ICR_genes = list(c(ICR_genes = ICR_genes))

# Analysis
Gene.set = get(Gene_set)
Expression.matrix = filtered.norm.RNAseqData
all_genes_in_set = unname(unlist(Gene.set))

available_genes = intersect(all_genes_in_set, rownames(Expression.matrix))
unavailabe_genes = all_genes_in_set[-which(all_genes_in_set %in% available_genes)]

ES = gsva(Expression.matrix, Gene.set, method="ssgsea")

# z score matrix
ESz = ES 
for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
}

# Annotation for plot
sam.info = data.frame(matrix(nrow = ncol(Expression.matrix), ncol = 0), row.names = colnames(Expression.matrix))
sam.info$ICR_cluster = table_cluster_assignment$ICR_HML[match(rownames(sam.info), rownames(table_cluster_assignment))]

od = hclust(dist(ESz))
nm = rownames(ESz)

col_fun = circlize::colorRamp2(c(min(ESz), 0, max(ESz)), c("blue", "white", "red"))

ha_column = HeatmapAnnotation(df = data.frame(ICR = sam.info$ICR_cluster),
                              show_annotation_name = TRUE,
                              col = list(ICR = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")),
                              show_legend = TRUE
)

#Heatmap
png(paste0("./Figures/Heatmaps/", Gene_set, "_Heatmaps/TCGA_COAD_Biolinks_", Gene_set, "_Heatmap.png"), res = 600, width = 13, height = 12, units = "in")

Heatmap(ESz, 
        name = "z scores", 
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_title_gp = gpar(fontsize = 0.1),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        top_annotation = ha_column,
        column_title = paste0(gsub("_", " ", Gene_set), " enrichment scores"),
        show_heatmap_legend = TRUE,
        row_names_max_width = unit(4, "in")
)

dev.off()

save(ES, file = paste0("./Analysis/Deconvolution_and_GSEA/", Gene_set, "_ES.Rdata"))
