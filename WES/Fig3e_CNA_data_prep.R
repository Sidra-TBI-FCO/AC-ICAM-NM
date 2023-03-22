
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

amp = read.table("./Processed_Data/WES/CNA/amp_genes.conf_95.txt", stringsAsFactors = FALSE,
                 fill = TRUE, sep = "\t")
del = read.table("./Processed_Data/WES/CNA/del_genes.conf_95.txt", stringsAsFactors = FALSE,
                 fill = TRUE, sep = "\t")
seg_file = read.table("./Processed_Data/WES/CNA/298SampleMerged.seg", stringsAsFactors = FALSE, 
                      sep = "\t", header = TRUE)

df_CNA = read.table("./Processed_Data/WES/CNA/all_thresholded.by_genes.txt", stringsAsFactors = FALSE,
                                      sep = "\t", header = TRUE)
load("./Overview Sample Stocks/Meta_data/Excluded_patients.Rdata")
load("./Overview Sample Stocks/Overview_available_data_per_patient/Old/Overview_21_September_2020_update.Rdata")

# Prepare data
WES_patients = Overview$Patient[which(Overview$WES == "yes" & Overview$Exclusion == "no")]

colnames(df_CNA)
colnames(df_CNA) = gsub("X", "", colnames(df_CNA))
colnames(df_CNA) = gsub("WES_COAD_LUMC_SIDRA_", "", colnames(df_CNA))
colnames(df_CNA)[4:ncol(df_CNA)] = str_pad(colnames(df_CNA)[4:ncol(df_CNA)], pad = "0", 4)

# Make vector of amplified and deleted genes
amp_genes = as.character(amp[5:nrow(amp), 2])

for (i in 3:ncol(amp)){
  amp_genes = c(amp_genes, as.character(amp[5:nrow(amp), i])) 
}
amp_genes = amp_genes[-which(amp_genes == "")]
amp_genes = amp_genes[-which(is.na(amp_genes))]

del_genes = as.character(del[5:nrow(del), 2])

for (i in 3:ncol(del)){
  del_genes = c(del_genes, as.character(del[5:nrow(del), i])) 
}
del_genes = del_genes[-which(del_genes == "")]
del_genes = del_genes[-which(is.na(del_genes))]

dir.create("./Analysis/WES/CNA", showWarnings = FALSE)
dir.create("./Analysis/WES/CNA/001_Data_prep", showWarnings = FALSE)

save(del_genes, amp_genes, file = "./Analysis/WES/CNA/001_Data_prep/del_amp_genes_vectors.Rdata")
save(df_CNA, file = "./Analysis/WES/CNA/001_Data_prep/CNA_All_thresholded_by_genes.Rdata")
save(seg_file, file = "./Analysis/WES/CNA/001_Data_prep/CNA_Seg_file.Rdata")

df_CNA = df_CNA[,c(1, 2, 3, which(substring(colnames(df_CNA), 1, 3) %in% WES_patients))]
gene_data = df_CNA[,1:3]
rownames(df_CNA) = df_CNA$Gene.Symbol
df_CNA = df_CNA[,-c(1:3)]

seg_file$Sample =  gsub("WES_COAD_LUMC_SIDRA_", "", seg_file$Sample)
seg_file$Sample = str_pad(seg_file$Sample, pad = "0", 4)
seg_file = seg_file[which(substring(seg_file$Sample, 1, 3) %in% WES_patients),]
length(unique(seg_file$Sample))

save(gene_data, df_CNA, file = "./Analysis/WES/CNA/001_Data_prep/CNA_281_patients_thresholded_by_genes.Rdata")
save(seg_file, file = "./Analysis/WES/CNA/001_Data_prep/CNA_281_patients_Seg_file.Rdata")
write.table(seg_file, file = "./Analysis/WES/CNA/001_Data_prep/CNA_281_patients_Seg_file.seg", sep = "\t", row.names = FALSE,
            quote = FALSE)

seg_file$Segment_Mean = log(seg_file$Segment_Mean, 2)
write.table(seg_file, file = "./Analysis/WES/CNA/001_Data_prep/CNA_281_patients_log2_Seg_file.seg", sep = "\t", row.names = FALSE,
            quote = FALSE)

Annotation = data.frame(Sample = unique(seg_file$Sample), CMS = NA, ICR = NA, MSI = NA,
                        POLE = NA, Hypermutation = NA)

Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed"
Annotation$CMS = Rfcms$RF.predictedCMS[match(Annotation$Sample, substring(rownames(Rfcms), 1, 4))]
Annotation$ICR = table_cluster_assignment$ICR_HML[match(Annotation$Sample, substring(rownames(table_cluster_assignment), 1, 4))]
Annotation$ICR[which(Annotation$ICR == "ICR Medium")] = "ICR Intermediate"
Annotation$MSI = MANTIS$MSI[match(Annotation$Sample, MANTIS$Sample_ID)]

# Add hypermutation status
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

Annotation$Hypermutation = frequency_df$Mutation_cat[match(substring(Annotation$Sample, 1, 3),
                                                                            frequency_df$Patient_ID)]

# Add POLE mutation
load("./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_clean_samples.Rdata")
length(unique(finalMafFiltered$Patient_ID))
POLE_mut = finalMafFiltered$Patient_ID[which(finalMafFiltered$Hugo_Symbol == "POLE")]

Annotation$POLE = "WT"
Annotation$POLE[which(substring(Annotation$Sample, 1, 3) %in% POLE_mut)] = "MUT"

Annotation = Annotation[, c("Sample", "Hypermutation", "POLE", "ICR", "CMS", "MSI")]

write.table(Annotation, file = "./Analysis/WES/CNA/001_Data_prep/New_Annotation_for_281_patients.txt", sep = "\t", row.names = FALSE,
            quote = FALSE)
