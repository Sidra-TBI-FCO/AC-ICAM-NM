
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ComplexHeatmap", "circlize")
ipak(required.packages)

# Set parameters
Immune_only = "Immune_only"

# Load data
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
sayaman_table = read.csv("./Processed_Data/External/Immune Traits Annotation/Sayaman_2021-mmc3.csv", stringsAsFactors = FALSE)

# Remove patients (from 366 -> 348 patients)
immune_sig_df = immune_sig_df[which(rownames(immune_sig_df) %in% colnames(RNASeq.QN.LOG2)),]

# For annotation table
sayaman_table$Sayaman.InternalLabel.Original = gsub("Sigs160.", "", sayaman_table$Sayaman.InternalLabel.Original)
sayaman_table$Sayaman.InternalLabel.Original = gsub("\\.", "_", sayaman_table$Sayaman.InternalLabel.Original)

annotation = data.frame(Full_label = colnames(immune_sig_df), Names = NA, Immune_trait_category = NA, Immune_module = NA)
annotation$Names = gsub(".*\\ - ", "", annotation$Full_label)
sayaman_table$Sayaman.FriendlyLabel[which(sayaman_table$Sayaman.FriendlyLabel == "Attractor Metagene - G HLA DPA1")] = "Attractor Metagene - G HLA-DPA1"
annotation$Immune_trait_category = sayaman_table$Sayaman.Annot.Figure.ImmuneCategory[match(annotation$Full_label,
                                                                                           sayaman_table$Sayaman.FriendlyLabel)]
annotation$Immune_module= sayaman_table$Sayaman.Annot.Figure.ImmuneModule[match(annotation$Full_label,
                                                                                sayaman_table$Sayaman.FriendlyLabel)]

annotation$Immune_trait_category[grep("ConsensusTME", annotation$Full_label)] = "ConsensusTME"
annotation$Immune_trait_category[which(annotation$Names %in% c("ISG RS", "IFNG GS"))] = "Expression Signature"
annotation$Immune_trait_category[which(is.na(annotation$Immune_trait_category))] = "Expression Signature"

if(Immune_only == "Immune_only"){
  colnames(immune_sig_df) == annotation$Full_label # check if all are in same order for next step
  colnames(immune_sig_df) = annotation$Names # annotation and immune-sig-df need to be in the same order
  annotation = annotation[which(annotation$Immune_trait_category %in% c("Attractor Metagene", "ConsensusTME", "Expression Signature")),]
  annotation$Full_label = annotation$Names
  immune_sig_df = immune_sig_df[, annotation$Full_label]
}

# Calculate correlation
mat_cor = cor(immune_sig_df, method = "pearson")
mat_cor = as.matrix(mat_cor)

dir.create("./Figures/Trimmed_p/048_Immune_traits_Pearson_correlation", showWarnings = FALSE)
#png(paste0("./Figures/Trimmed_p/048_Immune_traits_Pearson_correlation/048_2021_April_Pearson_correlation_heatmap_",
 #          Immune_only, ".png"), res = 600, 
  #  width = 13, height = 13, units = "in")
#hm = heatmap.2(mat_cor,
 #              # dendrogram control
  #             Rowv = TRUE,
   #            distfun = dist,
    #           hclustfun = hclust,
     #          dendrogram = c("both"),
      #         col = "bluered",
       #        trace = "none",
        #       margins=c(17,17),
         #      key = "T", density = "none")
# dev.off()

table(annotation$Immune_module)
#
column_ha = HeatmapAnnotation(Immune_trait_category = annotation$Immune_trait_category,
                              Immune_module = annotation$Immune_module,
                              show_annotation_name = FALSE,
                              simple_anno_size = unit("0.3", "cm"),
                              col = list(Immune_trait_category = c("Expression Signature" = "#6B9B78", "ConsensusTME" = "#BB2779",
                                                                   "Attractor Metagene" = "#556CA9"),
                                         Immune_module = c("IFN Response" = "#FFFEB0", "Lymphocyte Infiltration" = "#B5CADF",
                                                           "Macrophage/Monocyte" = "#C6DE93", "T-cell/Cytotoxic" = "#E8C083",
                                                           "TGF-b Response" = "#BDB0D1", "Unassigned" = "lightgrey", 
                                                           "Wound Healing" = "#ff6961")),
                              show_legend = FALSE)


col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "#ff2626"))

png(paste0("./Figures/Trimmed_p/048_Immune_traits_Pearson_correlation/048_April_2021_Pearson_correlation_Complex_heatmap_", Immune_only,".png"), 
    res = 600, width = 8, height = 8, units = "in")
HM = Heatmap(mat_cor,
             row_title_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5),
             column_title_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
             cluster_rows = TRUE, cluster_columns = T,
             show_column_names = T, 
             bottom_annotation = column_ha, 
             col = col_fun,
             show_heatmap_legend = TRUE,
             #col = color,
             # heatmap_legend_param =list(title_gp=gpar(fontsize=10, fontface="bold"),legend_width=unit(8,"cm"),legend_position = "left"),
             #row_names_max_width = unit(5, "cm")
             # theme(legend.position = "none")
             row_names_max_width = unit(6, "in")
)
HM = draw(HM, heatmap_legend_side = "left")
dev.off()

pdf(paste0("./Figures/Trimmed_p/048_Immune_traits_Pearson_correlation/048_April_2021_Pearson_correlation_Complex_heatmap_", Immune_only,".pdf"), 
    width = 8, height = 8)
HM = draw(HM, heatmap_legend_side = "left")
dev.off()

order = rownames(mat_cor)[row_order(HM)]

df = data.frame(Full_name = order, Module = NA)
df$Module[1:29] = "M1"
df$Module[30:57] = "M2"
df$Module[58:72] = "M3"
df$Module[73:82] = "M4"
df$Module[83:95] = "M5"
df$Module[96:100] = "M6"
df$Module[101:103] = "M7"

annotation$Module_Roelands = df$Module[match(annotation$Full_label, df$Full_name)]

dir.create("./Analysis/Trimmed_p/048_Correlation_plot_order_pathways", showWarnings = FALSE)
save(order, file = "./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/048_2021_Correlation_plot_order_pathways.Rdata")

save(df, annotation, file = "./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/Immune_signatures_clusters_Roelands.Rdata")

load("./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/Immune_signatures_clusters_Roelands.Rdata")
table(df$Module)
i = 1

results = data.frame(Modules = unique(df$Module), Average_correlation = NA, sd = NA)

i=7
for (i in 1:7){
  Module = paste0("M", i, sep = "")
  test = mat_cor[df$Full_name[which(df$Module == Module)],df$Full_name[which(df$Module == Module)]]
  results$Average_correlation[which(results$Modules == Module)] = mean(test)
  results$sd[which(results$Modules == Module)] = sd(test)
}

write.csv(results, file = "./Analysis/Trimmed_p/048_Correlation_plot_order_pathways/048_Mean_correlation_per_immune_module.csv",
          row.names = FALSE)



