
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr", "ggplot2")
ipak(required.packages)

# Load data
Sidra_LUMC = read.csv("./Analysis/WES/011_Compare_with_TCGA/Oct_2021_Sidra_LUMC_frequency_compared_with_TCGA_NHS_PFHS.csv",
                      stringsAsFactors = FALSE)

# Only keep genes that are likely involved in cancer development
# QGP genes collection (n=1226)
# No artifact genes
Filtered = Sidra_LUMC[which(Sidra_LUMC$Gene.in.QGP.collection == "yes"),]
Filtered = Filtered[which(Filtered$Artifact.gene == "no"),]

# Select those that have not been defined as colon oncogenic mediator by Colaprico et al
# Select those genes that have not been found in Bailey COADREAD
Filtered_new = Filtered
Filtered_new = Filtered_new[which(Filtered_new$Gene.in.Colaprico.COAD == "no"),]
Filtered_new = Filtered_new[which(Filtered_new$Gene.in.Bailey.COADREAD == "no"),]

# Filter out genes that were significant in MutSig analysis of Grasso et al
#Filtered_new = Filtered_new[-which(Filtered_new$Grasso_MutSig_in_MSI_H_group == "yes"),]
#Filtered_new = Filtered_new[-which(Filtered_new$Grasso_MutSig_in_MSS_group == "yes"),]

#  Filter out genes that were significant in MutSig analysis of Giannakis et al
#Filtered_new = Filtered_new[which(is.na(Filtered_new$Giannakis_MutSig_Percentage_in_NHS_HPFS_nonhypermutated)),]
#Filtered_new = Filtered_new[which(is.na(Filtered_new$Giannakis_MutSig_Percentage_in_NHS_HPFS_hypermutated)),]

# Keep genes that have a frequency of <5% in NHS HPFS Colon Cancer cohort
Filtered_new = Filtered_new[which(Filtered_new$Giannakis..Percentage.in.NHS.HPFS.Colon.Cancer < 5),]

# Keep genes that have a frequency of <5% in TCGA
Filtered_new = Filtered_new[which(Filtered_new$Percentage.in.TCGA.COAD < 5),]

# Keep genes that have a frequency of <5% in COSMIC
#Filtered_new = Filtered_new[which(Filtered_new$Percentage_in_COSMIC_Colon_Adenocarcinoma < 5),]

# Keep genes that have a frequency >5% in Sidra-LUMC
Filtered_new = Filtered_new[which(Filtered_new$Percentage.in.Sidra.LUMC.cohort > 5),]

# Keep genes that are at least twice more increased in Sidra-LUMC compared to TCGA-COAD
Filtered_new$fold_increase_in_Sidra_LUMC = Filtered_new$Percentage.in.Sidra.LUMC.cohort / Filtered_new$Percentage.in.TCGA.COAD
Filtered_new = Filtered_new[which(Filtered_new$fold_increase_in_Sidra_LUMC > 2),]

# Annotate Filtered with a column that specifies if a gene is new
Filtered$New = "Old"
Filtered$New[which(Filtered$Hugo.Symbol %in% Filtered_new$Hugo.Symbol)] = "New"

plot_df = Filtered[which(Filtered$Percentage.in.Sidra.LUMC.cohort > 5),]
plot_df = plot_df[-c(1:21),]
plot_df = plot_df[-which(plot_df$Hugo.Symbol == "POLE"),]
plot = ggplot(plot_df, aes(x = Hugo.Symbol, y = Percentage.in.Sidra.LUMC.cohort,
                           fill = New)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("New" = "#FF6273", "Old" = "#084BDD")) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 12.5,
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Mutation \n frequency (%)") +
  xlab("") +
  scale_x_discrete(limits= plot_df$Hugo.Symbol) +
  ylim(0, 15)

pdf("./Figures/WES/011_Frequently_mutated_genes/Fig3c__Barplot_frequently_mut_genes.pdf",
    width = 32, height = 1.8)
plot(plot)
dev.off()
