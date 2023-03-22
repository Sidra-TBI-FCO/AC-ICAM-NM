#################################################################
###
### Create RNASEQ based Deconvolution
### using Bindea's signatures and ssGSEA
###
#################################################################

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("GSVA","ComplexHeatmap", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
Gene_set = "Selected.pathways" #"Selected.pathways.Sidra.LUMC" # "Cholesterol_list" "MSigDB_list" "Angelove_ORIG" "Bindea_ORIG" "ICR" "Bindea_PATRICK" "Thorsson" "Prat_ORIG_NANO" "Trajanoski_ORIG"  "Thorsson" "Benci_list" "Selected.pathways" "TIS" "Specific genes"
subset = ""

# Create folders and log file
dir.create("./Figures",showWarnings = FALSE)
dir.create("./Figures/Trimmed_p/Heatmaps", showWarnings = FALSE)
dir.create(paste0("./Figures/Trimmed_p/Heatmaps/", Gene_set, "_Heatmaps"), showWarnings = FALSE)
dir.create("./Analysis/Trimmed_p/Deconvolution_and_GSEA", showWarnings = FALSE)

#load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load(paste0(toolbox.path, "/GSEA list/immune.gene.lists.v4.Rdata"))
load(paste0(toolbox.path, "/Selected Pathways/Selected.pathways.3.4.RData"))
load(paste0(toolbox.path, "/GSEA list/GSEA_Module_list_Gen3.Rdata"))
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))
load(paste0(toolbox.path, "/GSEA list/DDR.gene.list_v2.Rdata"))
load(paste0(toolbox.path, "/TIS nanostring genes/18_gene_TIS.Rdata"))
load(paste0(toolbox.path, "/GSEA list/MSigDB_list.Rdata"))
load(paste0(toolbox.path, "/GSEA list/Cholesterol_list.Rdata"))
load("./Sidra_LUMC_specific_R_tools/GSEA_list/Immune_gene_signatures_v1.Rdata")
load("./Sidra_LUMC_specific_R_tools/GSEA_list/Selected_pathways_v1.Rdata")
#BRCAness = read.csv("./Processed_Data/External/Gene_collections/Konstantinopoulus_BRCAness_Signature_Supplementary_Table_1.csv", stringsAsFactors = FALSE)
#BRCAness_genes = BRCAness$Gene.Symbol

ICR <-list()
list.item = list(ICR_genes)
names(list.item) = "ICR_genes"
ICR = c(ICR, list.item)

TIS <-list()
list.item = list(TIS_genes)
names(list.item) = "TIS_genes"
TIS = c(TIS, list.item)
#BRCAness = list()
#list.item = list(BRCAness_genes)
#names(list.item) = "BRCAness"
#BRCAness = c(BRCAness, list.item)

if(Gene_set == "single_genes_list"){
  single_genes_list = list()
  list.item = list(c("TGFB1"))
  names(list.item) = "TGFB1"
  list.item2 = list(c("IFNG"))
  names(list.item2) = "IFNG"
  single_genes_list = c(single_genes_list, list.item, list.item2)
}

if(subset == ""){
  load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.clean.dataset.EDAseq.QN.HPC.Rdata")
}
if(subset == "including_all"){
  load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.Complete.dataset.EDAseq.QN.HPC.Rdata")
}

# Analysis
Gene.set = get(Gene_set)
#Gene.set = c(Gene.set, list.item)
Expression.matrix = RNASeq.QN.LOG2
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
png(paste0("./Figures/Trimmed_p/Heatmaps/", Gene_set, "_Heatmaps/ssGSEA_JSREP_",subset, Gene_set, "_Heatmap.png"), 
    res = 600, width = 9, height = 6, units = "in")

Heatmap(ESz, 
        name = "z scores", 
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_title_gp = gpar(fontsize = 0.1),
        column_names_gp = gpar(fontsize = 0),
        row_names_gp = gpar(fontsize = 10),
        top_annotation = ha_column,
        column_title = paste0(gsub("_", " ", Gene_set), " enrichment scores"),
        show_heatmap_legend = TRUE,
        row_names_max_width = unit(4, "in")
)

dev.off()

save(ES, ESz, file = paste0("./Analysis/Trimmed_p/Deconvolution_and_GSEA/",subset, Gene_set, "_ES.Rdata"))
