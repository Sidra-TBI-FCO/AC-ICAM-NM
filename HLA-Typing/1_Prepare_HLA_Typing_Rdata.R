
## Prepare HLA-typing merged file

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr")
ipak(required.packages)

# # Prepare RNASeq HLA Typing result
Optitype_RNA = read.csv("./Processed_Data/HLA-Typing/RNA_HLA_Optitype/merged_all.csv", stringsAsFactors = FALSE)
Optitype_RNA = Optitype_RNA[-which(Optitype_RNA$sample == "sample"),]
Optitype_RNA$index = NULL
Optitype_RNA$score = NULL
colnames(Optitype_RNA) = c("sample", "locus", "alleles_RNASeq")
Optitype_RNA$Patient_ID = substring(Optitype_RNA$sample, 1, 3)
Optitype_RNA$Patient_ID_Locus = paste(Optitype_RNA$Patient_ID, Optitype_RNA$locus, sep = "_")
Optitype_RNA$Tissue_type = substring(Optitype_RNA$sample, 4, 4)
Optitype_RNA$Tissue_type[which(Optitype_RNA$Tissue_type == "L")] = "LM"
table(Optitype_RNA$Tissue_type, exclude = NULL)

# Split according to tissue type (seperate column in final results table)
Optitype_RNA_T = Optitype_RNA[which(Optitype_RNA$Tissue_type == "T"),]
Optitype_RNA_LM = Optitype_RNA[which(Optitype_RNA$Tissue_type == "LM"),]

# Prepare WES HLA Typing result
Optitype_WES = read.csv("./Processed_Data/HLA-typing/WES_HLA_Optitype/Final_HLA_Optitype_results_WES.csv", stringsAsFactors = FALSE)
Optitype_WES$expected = NULL
Optitype_WES$validates = NULL
colnames(Optitype_WES) = c("sample", "locus", "alleles_WES")
Optitype_WES$sample = gsub("WES_COAD_LUMC_SIDRA_", "", Optitype_WES$sample)
Optitype_WES$sample = substring(Optitype_WES$sample, 1, nchar(Optitype_WES$sample) - 10)
Optitype_WES$sample = str_pad(Optitype_WES$sample, 4, pad = "0")
Optitype_WES$Patient_ID = substring(Optitype_WES$sample, 1, 3)
Optitype_WES$Tissue_type = substring(Optitype_WES$sample, 4, 4)
Optitype_WES$Tissue_type[which(Optitype_WES$Tissue_type == "L")] = "LM"
Optitype_WES$Patient_ID_Locus = paste(Optitype_WES$Patient_ID, Optitype_WES$locus, sep = "_")
table(Optitype_WES$Tissue_type, exclude = NULL)

# Split according to tissue type (seperate column in final results table)
Optitype_WES_T = Optitype_WES[which(Optitype_WES$Tissue_type == "T"),]
Optitype_WES_N = Optitype_WES[which(Optitype_WES$Tissue_type == "N"),]
Optitype_WES_LM = Optitype_WES[which(Optitype_WES$Tissue_type == "LM"),]

# Combine the data, 3 single rows per patient (for HLA-A-B-C)
all_patients = unique(c(Optitype_RNA$Patient_ID, Optitype_WES$Patient_ID))
all_results_df = data.frame(Patient_ID = rep(all_patients,3),
                            Locus = c(rep("A",length(all_patients)), rep("B", length(all_patients)), rep("C", length(all_patients))),
                            Patient_ID_Locus = NA,
                            Alleles_T_RNASeq = NA, Alleles_T_WES = NA, 
                            Alleles_N_WES = NA, Alleles_LM_RNASeq = NA,
                            Alleles_LM_WES = NA)
all_results_df$Patient_ID_Locus = paste(all_results_df$Patient_ID, all_results_df$Locus, sep = "_")

# Add RNASeq data (T and LM in seperate column)
all_results_df$Alleles_T_RNASeq = Optitype_RNA_T$alleles_RNASeq[match(all_results_df$Patient_ID_Locus,
                                                                      Optitype_RNA_T$Patient_ID_Locus)]

all_results_df$Alleles_LM_RNASeq = Optitype_RNA_LM$alleles_RNASeq[match(all_results_df$Patient_ID_Locus,
                                                                        Optitype_RNA_LM$Patient_ID_Locus)]
                                                                    

# Add WES data (T, N, and LM in seperate columns)
all_results_df$Alleles_T_WES = Optitype_WES_T$alleles_WES[match(all_results_df$Patient_ID_Locus,
                                                                Optitype_WES_T$Patient_ID_Locus)]


all_results_df$Alleles_N_WES = Optitype_WES_N$alleles_WES[match(all_results_df$Patient_ID_Locus,
                                                                Optitype_WES_N$Patient_ID_Locus)]


all_results_df$Alleles_LM_WES = Optitype_WES_LM$alleles_WES[match(all_results_df$Patient_ID_Locus,
                                                                  Optitype_WES_LM$Patient_ID_Locus)]


# Order elements in each of the columns
tmp = lapply(str_split(all_results_df$Alleles_T_RNASeq, ";"), sort)
new = lapply(tmp, paste, collapse = ";")
all_results_df$Alleles_T_RNASeq = as.character(new)

tmp = lapply(str_split(all_results_df$Alleles_T_WES, ";"), sort)
new = lapply(tmp, paste, collapse = ";")
all_results_df$Alleles_T_WES = as.character(new)

tmp = lapply(str_split(all_results_df$Alleles_N_WES, ";"), sort)
new = lapply(tmp, paste, collapse = ";")
all_results_df$Alleles_N_WES = as.character(new)

tmp = lapply(str_split(all_results_df$Alleles_LM_WES, ";"), sort)
new = lapply(tmp, paste, collapse = ";")
all_results_df$Alleles_LM_WES = as.character(new)

tmp = lapply(str_split(all_results_df$Alleles_LM_RNASeq, ";"), sort)
new = lapply(tmp, paste, collapse = ";")
all_results_df$Alleles_LM_RNASeq = as.character(new)

all_results_df[ all_results_df==""] = NA

dir.create("./Analysis/HLA Typing/Optitype_merged_result", showWarnings = FALSE)
save(all_results_df, file = "./Analysis/HLA Typing/Optitype_merged_result/Results_Optitype_RNASeq_WES_by_Patient_ID.Rdata")
write.csv(all_results_df, file = "./Analysis/HLA Typing/Optitype_merged_result/Results_Optitype_RNASeq_WES_by_Patient_ID.csv")
