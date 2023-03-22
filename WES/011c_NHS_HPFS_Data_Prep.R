
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("dplyr")
ipak(required.packages)

# Load data
anatomic_site = read.csv("./Processed_Data/External/TableS1_Giannakis_et_al_supplementary_tables.csv", 
                         stringsAsFactors = FALSE)

MAF_long = read.csv("./Processed_Data/External/TableS3_Giannakis_et_al_supplementary_tables.csv",
                    stringsAsFactors = FALSE)

# Checks
patients_anatomic = anatomic_site$Individual.ID
MAF_long$Patient_ID = gsub("T", "", MAF_long$Tumor_Sample_Barcode)
patients = unique(MAF_long$Patient_ID)
patients %in% patients_anatomic # Should all be TRUE
patients_anatomic %in% patients # Should all be TRUE

table(MAF_long$Variant_Classification)
MAF_non_silent = MAF_long[which(MAF_long$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", 
                                                                       "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")),]
table(MAF_non_silent$Variant_Classification)

colon_patients = anatomic_site$Individual.ID[-which(anatomic_site$Tumor.Site == "rectum")]
MAF_non_silent = MAF_non_silent[which(MAF_non_silent$Patient_ID %in% colon_patients),] 

genes = unique(MAF_non_silent$Hugo_Symbol)
patients = unique(MAF_non_silent$Patient_ID)

matrix = matrix(nrow = length(genes), ncol = length(patients))
rownames(matrix) = genes
colnames(matrix) = patients

df = data.frame(Gene = genes, n_mut = NA)
i=1
for (i in 1:length(genes)){
  gene = genes[i]
  MAF_sub = MAF_non_silent[which(MAF_non_silent$Hugo_Symbol == gene),]
  patients_mut = unique(MAF_sub$Patient_ID)
  df$n_mut[which(df$Gene == gene)] = length(patients_mut)
}

df$Frequency = df$n_mut / length(patients)

save(df, file = "./Analysis/WES/011_Compare_with_TCGA/NHS_HPFS_Data_Frequency_Nonsilent_mutations_in_Colon_only_482_patients.Rdata")
