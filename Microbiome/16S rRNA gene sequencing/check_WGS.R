# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# load packages
ipak(c("dplyr", "reshape2", "ggplot2", "VennDiagram", "openxlsx", "stringr", "venn"))

# load data
load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/Microbiome/Dohlman_TCGA_COAD_no_duplicate_filtered_per_10_RA_0.01.Rdata")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_N_246_samples_based_on_normal.Rdata")
normal = Genus_full_abundance
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
fil.WGS = read.csv("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/WGS_validation/Genera_filtered_Relative_abundance_matrix_WGS_167.csv")
tcga.phylum = read.table("../NGS_Data_TCGA_COAD_Jessica/Microbiome data TCMA/bacteria.sample.relabund.phylum.txt", sep = '\t', header = T)
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")

#rownames(Genus_full_abundance)
Genus_full_abundance = as.data.frame(Genus_full_abundance)
Genus_full_abundance$genus = gsub(".*\\D_5__", "", rownames(Genus_full_abundance))

tumor.genera = as.data.frame(rownames(Genus_full_abundance))
Genus_full_abundance$genus[which(Genus_full_abundance$genus == "Ruminococcus 2")] = "Ruminococcus"
Genus_full_abundance$genus[which(rownames(Genus_full_abundance) == "D_0__Bacteria D_1__Bacteroidetes D_2__Bacteroidia D_3__Bacteroidales D_4__Prevotellaceae D_5__uncultured")] = "Prevotellaceae_uncultured"
Genus_full_abundance$genus[which(rownames(Genus_full_abundance) == "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Lachnospiraceae D_5__uncultured")] = "Lachnospiraceae_uncultured"
Genus_full_abundance$genus[which(rownames(Genus_full_abundance) == "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__uncultured")] = "Ruminococcaceae_uncultured"
Genus_full_abundance$genus[which(rownames(Genus_full_abundance) == "D_0__Bacteria D_1__Firmicutes D_2__Erysipelotrichia D_3__Erysipelotrichales D_4__Erysipelotrichaceae D_5__uncultured")] = "Erysipelotrichaceae_uncultured"
Genus_full_abundance$genus[which(rownames(Genus_full_abundance) == "D_0__Bacteria D_1__Proteobacteria D_2__Deltaproteobacteria D_3__Desulfovibrionales D_4__Desulfovibrionaceae D_5__uncultured")] = "Desulfovibrionaceae_uncultured"

# 16s tumor phylum
Genus_full_abundance$phylum = gsub("\\D_2__.*", "", rownames(Genus_full_abundance))
Genus_full_abundance$phylum = gsub(".*\\D_1__", "", Genus_full_abundance$phylum)

# TCGA
micro = tcga_bacteria_df

# phylum
rownames(tcga.phylum) = tcga.phylum$name
tcga.phylum$name = NULL

colnames(tcga.phylum) = gsub("\\.", "-", colnames(tcga.phylum))
tcga.phylum = tcga.phylum[, colnames(micro)]

### WGS
WGS = Genus_WGS
#rownames(WGS) = WGS$genus
rownames(WGS) = gsub("g__", "", rownames(WGS))  # genus
#rownames(WGS)[which(rownames(WGS) == "Clostridiales_Family_XIII_Incertae_Sedis_unclassified")] = "Clostridiales Family XIII Incertae Sedis unclassified"

#WGS$genus = NULL
#WGS$All_taxa = NULL

rownames(Phylum_WGS) = gsub("p__", "", rownames(Phylum_WGS))  # phylum

# normal
normal = as.data.frame(normal)
normal$genus = gsub(".*\\D_5__", "", rownames(normal))
normal$genus[which(normal$genus == "Ruminococcus 2")] = "Ruminococcus"
normal$genus[which(rownames(normal) == "D_0__Bacteria D_1__Bacteroidetes D_2__Bacteroidia D_3__Bacteroidales D_4__Prevotellaceae D_5__uncultured")] = "Prevotellaceae_uncultured"
normal$genus[which(rownames(normal) == "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Lachnospiraceae D_5__uncultured")] = "Lachnospiraceae_uncultured"
normal$genus[which(rownames(normal) == "D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__uncultured")] = "Ruminococcaceae_uncultured"
normal$genus[which(rownames(normal) == "D_0__Bacteria D_1__Firmicutes D_2__Erysipelotrichia D_3__Erysipelotrichales D_4__Erysipelotrichaceae D_5__uncultured")] = "Erysipelotrichaceae_uncultured"
normal$genus[which(rownames(normal) == "D_0__Bacteria D_1__Proteobacteria D_2__Deltaproteobacteria D_3__Desulfovibrionales D_4__Desulfovibrionaceae D_5__uncultured")] = "Desulfovibrionaceae_uncultured"

# phylum
normal$phylum = gsub("\\D_2__.*", "", rownames(normal))
normal$phylum = gsub(".*\\D_1__", "", normal$phylum)

# remove space
Genus_full_abundance$phylum = str_trim(Genus_full_abundance$phylum, "both")
normal$phylum = str_trim(normal$phylum, "both")

# Venn genus unfil WGS

venn.plot = 
  list('16S tumor' = Genus_full_abundance$genus, '16S normal' = normal$genus ,'WGS AC-ICAM' = rownames(WGS), TCGA = rownames(micro))

#svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/venndiagram/all_filtered_genera_overlapp_16S_tumor_normal_WGS_TCGA.svg"), width = 3, height = 3)

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#fa9b91","#78acff" ,"#A0FFA0", "gold"), ilcs = 1, sncs = 0.7);

dev.off()

common = Reduce(intersect,list(Genus_full_abundance$genus, normal$genus ,rownames(Genus_WGS), rownames(micro)))

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/venndiagram/genera_overlapp_16S_tumor_normal_WGS_TCGA.png"), width = 3, height = 3, res = 600, units = "in")

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#fa9b91","#78acff" ,"#A0FFA0", "gold"), ilcs = 1, sncs = 0.7);

dev.off()

################

# Venn phylum 

venn.plot = 
  list('16S tumor' = Genus_full_abundance$phylum, '16S normal' = normal$phylum,'WGS AC-ICAM' = rownames(Phylum_WGS), TCGA = rownames(tcga.phylum))

svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/venndiagram/phylum_overlapp_16S_tumor_normal_WGS_TCGA.svg"), width = 3, height = 3)

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#fa9b91","#78acff" ,"#A0FFA0", "gold"), ilcs = 1, sncs = 0.7);

dev.off()

Reduce(intersect,list(Genus_full_abundance$phylum, normal$phylum, rownames(Phylum_WGS), tcga.phylum$name))

png(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/venndiagram/phylum_overlapp_16S_tumor_normal_WGS_TCGA.png"), width = 3, height = 3, res = 600, units = "in")

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#fa9b91","#78acff" ,"#A0FFA0", "gold"), ilcs = 1, sncs = 0.7);

dev.off()

# phylum.tcga = tcga.phylum
# rownames(phylum.tcga) = phylum.tcga$name
# phylum.tcga$name = NULL
# 
# phylum.tcga = as.data.frame(t(phylum.tcga))




# WGS$phylum = gsub(".*\\|p__", "", WGS$All_taxa)
# WGS$phylum = gsub("\\|c__.*", "", WGS$phylum) 
# 
# table(WGS$phylum)

common = Reduce(intersect, list(Genus_full_abundance$genus, rownames(WGS))) ## unfiltered 

fil.WGS$genus = gsub("g__", "", fil.WGS$genus)
fi.comm = Reduce(intersect, list(Genus_full_abundance$genus, fil.WGS$genus)) ## filtered 

# Venn genus fil WGS

venn.plot = 
  list('16S tumor' = Genus_full_abundance$genus, '16S normal' = normal$genus ,'WGS AC-ICAM' = fil.WGS$genus, TCGA = rownames(micro))

#svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/venndiagram/all_filtered_genera_overlapp_16S_tumor_normal_WGS_TCGA.svg"), width = 3, height = 3)

venn.result =
  venn(venn.plot, ilabels = TRUE, 
       zcolor = c("#fa9b91","#78acff" ,"#A0FFA0", "gold"), ilcs = 1, sncs = 0.7);

dev.off()

common = Reduce(intersect,list(Genus_full_abundance$genus, normal$genus ,rownames(Genus_WGS), rownames(micro)))

##############################################################################

check.table = read.xlsx("C:/Users/eahmed2/Desktop/Copy of Supplementary_Tables_2-12.xlsx")

common1 = Reduce(intersect,list(fil.WGS$genus, check.table$AC.ICAM.WGS.Genus))

check.table$AC.ICAM.WGS.Genus[-which(check.table$AC.ICAM.WGS.Genus %in% fil.WGS$genus)]


table.list = as.data.frame(check.table$AC.ICAM.WGS.Genus)
fil.list = as.data.frame(fil.WGS$genus)


common.fil = Reduce(intersect, list(Genus_full_abundance$genus, fil.WGS$genus))
common.unfil = Reduce(intersect, list(Genus_full_abundance$genus, rownames(WGS)))

