
## Annotate results with clinical data and ICR

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2")
ipak(required.packages)

# Load data
TCR_Overview = read.csv("./Processed_Data/TCR/TCR SampleOverview/SampleOverview_kit_1_2_3.csv", stringsAsFactors = FALSE)
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
ICR_primaries = table_cluster_assignment
load("./Analysis/Trimmed_p/ICR Consensus Clustering/With_metastasis/JSREP_ICR_cluster_assignment_k2-6.Rdata")
ICR_with_meta = table_cluster_assignment
rm(list = "table_cluster_assignment")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

# Annotation
rownames(ICR_with_meta) = gsub("_P", "", rownames(ICR_with_meta))
rownames(ICR_primaries) = gsub("_P", "", rownames(ICR_primaries))

TCR_Overview$X = NULL
TCR_Overview$Patient_ID = substring(TCR_Overview$sample_name, 1, 3)
TCR_Overview$Tissue = substring(TCR_Overview$sample_name, 4, nchar(TCR_Overview$sample_name))
TCR_Overview$Tissue[which(TCR_Overview$Tissue %in% c("LM1", "LM2"))] = "Liver metastasis"
TCR_Overview$Tissue[which(TCR_Overview$Tissue %in% c("T"))] = "Primary tumor"
TCR_Overview$Tissue[which(TCR_Overview$Tissue %in% c("N"))] = "Normal colon"
TCR_Overview$Tissue = factor(TCR_Overview$Tissue, levels = c("Normal colon", "Primary tumor", "Liver metastasis"))

TCR_Overview$ICRscore = ICR_with_meta$ICRscore[match(TCR_Overview$sample_name, rownames(ICR_with_meta))]
TCR_Overview$ICRscore[which(TCR_Overview$sample_name == "030T")] = ICR_primaries$ICRscore[which(substring(rownames(ICR_primaries), 1, 4) == "030T")]

TCR_Overview$HLM_cluster = ICR_primaries$ICR_HML[match(TCR_Overview$sample_name, rownames(ICR_primaries))]
TCR_Overview$HLM_cluster_meta = ICR_with_meta$ICR_HML[match(TCR_Overview$sample_name, rownames(ICR_with_meta))]
TCR_Overview$HLM_cluster = factor(TCR_Overview$HLM_cluster, levels = c("ICR High", "ICR Medium", "ICR Low"))
TCR_Overview$HLM_cluster_meta = factor(TCR_Overview$HLM_cluster_meta, levels = c("ICR High", "ICR Medium", "ICR Low"))

TCR_Overview$OS.Time = clinical_data$OS.Time[match(TCR_Overview$Patient_ID,clinical_data$Patient_ID)]
TCR_Overview$OS.Status = clinical_data$OS.Status[match(TCR_Overview$Patient_ID,clinical_data$Patient_ID)]
TCR_Overview$primary_tumour_anatomic_site = clinical_data$tumour_anatomic_site[match(TCR_Overview$Patient_ID,clinical_data$Patient_ID)]
TCR_Overview$primary_tumour_anatomic_site = factor(TCR_Overview$primary_tumour_anatomic_site, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                                                         "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                                                         "rectosigmoideum"))
TCR_Overview$primary_tumor_side = NA
TCR_Overview$primary_tumor_side[which(TCR_Overview$primary_tumour_anatomic_site %in% c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                               "colon transversum"))] = "Right sided"
TCR_Overview$primary_tumor_side[which(TCR_Overview$primary_tumour_anatomic_site %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                               "rectosigmoideum"))] = "Left sided"
TCR_Overview$primary_tumor_side = factor(TCR_Overview$primary_tumor_side, levels = c("Right sided", "Left sided"))
TCR_Overview$primary_tumour_ajcc_pathologic_stage = clinical_data$ajcc_pathologic_tumor_stage[match(TCR_Overview$Patient_ID,clinical_data$Patient_ID)]
TCR_Overview$fraction_unique_clonotypes = TCR_Overview$productive_rearrangements / TCR_Overview$productive_templates

TCR_Overview = TCR_Overview[-which(TCR_Overview$Patient_ID %in% excluded_df$Patient_ID),]

dir.create("./Analysis/TCR", showWarnings = FALSE)
dir.create("./Analysis/TCR/0_1_Annotated_data", showWarnings = FALSE)
save(TCR_Overview, file = "./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
write.csv(TCR_Overview, file = "./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.csv", row.names = FALSE)
