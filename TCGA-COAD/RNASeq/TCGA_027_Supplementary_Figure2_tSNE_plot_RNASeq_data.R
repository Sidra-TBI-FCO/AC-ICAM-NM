
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("ggplot2", "Rtsne")                                                                   
ibiopak(required.bioconductor.packages)

# Set parameters
Group = "ICR_HML"

#load data
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")

#Analysis
RNASeq.QN.LOG2 = t(filtered.norm.RNAseqData)

set.seed(7)
#tsne_model_1 = Rtsne(RNASeq.QN.LOG2, check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)

#dir.create("./Analysis/027_tSNE_plot_RNASeq", showWarnings = FALSE)
#save(tsne_model_1, file = "./Analysis/027_tSNE_plot_RNASeq/RNASeq_tsne_model_3.Rdata")
load("./Analysis/027_tSNE_plot_RNASeq/RNASeq_tsne_model_3.Rdata")

df = table_cluster_assignment
df$Patient_ID = substring(rownames(df), 1, 12)
#df$tumor_morphology = clinical_data$Tumor_morphology[match(df$Patient_ID, clinical_data$Patient_ID)]
df$CMS = Rfcms$RF.predictedCMS[match(df$Patient_ID, substring(rownames(Rfcms), 1, 12))]
df$CMS[which(is.na(df$CMS))] = "Not predicted"
#df$tumor_anatomic_location = clinical_data$tumour_anatomic_site[match(df$Patient_ID, clinical_data$Patient_ID)]
#df$tumor_anatomic_location = factor(df$tumor_anatomic_location,
       #                             levels = c("ceceum", "colon ascendens",
        #                                       "flexura hepatica", "colon transversum",
         #                                      "flexura lienalis", "colon descendens",
         #                                      "colon sigmoideum", "rectosigmoideum"))
df = df[rownames(RNASeq.QN.LOG2),]

if(Group == "tumor_morphology"){
  colors = c("grey", "green", "blue", "purple", "red", "orange", "#FF3399",
             "darkgreen", "cyan", "darkblue", "black", "#009999", "darkred")
  names(colors) = unique(df$tumor_morphology)
}
if(Group == "CMS"){
  colors = c("CMS1" = "#FF9F21", "CMS2" = "#0074AF", "CMS3" = "#E97AA8", "CMS4" = "#009E74",
             "Not predicted" = "lightgrey")
}
if(Group == "ICR_HML"){
  colors = c("ICR High" = "red", "ICR Medium" = "green", "ICR Low" = "blue")
}
if(Group == "tumor_anatomic_location"){
  colors = c("green", "blue", "purple", "red", "orange", "#FF3399",
             "darkgreen", "red")
}


#dev.new()
#plot(tsne_model_1$Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
#text(tsne_model_1$Y, labels=df$tumor_morphology, col=colors[df$tumor_morphology])

d_tsne_1 = as.data.frame(tsne_model_1$Y) 

dir.create("./Figures/027_tSNE_plot_RNASeq", showWarnings = FALSE)
png(paste0("./Figures/027_tSNE_plot_RNASeq/v3_tSNE_by_", Group, ".png"), height = 3.8, width = 4.2,
    units = "in", res = 600)
plot= ggplot(d_tsne_1, aes(x=V1, y=V2,color=df[, Group])) +
  geom_point(aes(x=V1, y=V2,fill=df[, Group]),size=1.3) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  #scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=colors)+
  xlab("tSNE dimension 1") + ylab("tSNE dimension 2") +
  #ggtitle("Perplexity=15") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none"
  )
plot(plot)
dev.off()
