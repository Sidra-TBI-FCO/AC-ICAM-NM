
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "dplyr"))

# Load data
load("./Overview Sample Stocks/Overview_available_data_per_patient/Overview_26_April_2022.Rdata")
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Clean_TRB_clonality_341_patients.Rdata")
load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/WES/101_AC-ICAM_Calculate_GIE_Clean/frequency_df_GIE_AC_ICAM.Rdata")
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")
load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")

# Prepare data
df = Overview
colnames(df)[1] = "Patient_ID"
colnames(df)[5] = "Microbiome_Paired_TN"

# ICR
df$var01_ICRscore = table_cluster_assignment$ICRscore[match(df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
df$var02_ICR_cluster = table_cluster_assignment$ICR_HML[match(df$Patient, substring(rownames(table_cluster_assignment), 1, 3))]
df$var03_ICR_UP_DOWN = NA
median_ICR = median(df$var01_ICRscore)
df$var03_ICR_UP_DOWN[which(df$var01_ICRscore >=median_ICR)] = "UP"
df$var03_ICR_UP_DOWN[which(df$var01_ICRscore <median_ICR)] = "DOWN"
df$var03_ICR_UP_DOWN = factor(df$var03_ICR_UP_DOWN, levels = c("DOWN", "UP"))
# check
table(df$var02_ICR_cluster, df$var03_ICR_UP_DOWN, exclude = NULL)

# MiXCR TCR clonality
df$var04_MiXCR_clonality = MiXCR$MiXCR_Clonality[match(df$Patient_ID, substring(MiXCR$Sample_ID, 1, 3))]
df$var04_MiXCR_clonality[which(df$var04_MiXCR_clonality == -Inf)] = 0
df$var05_MiXCR_clonality_group = NA
median_MiXCR = median(MiXCR$MiXCR_Clonality, na.rm = TRUE)
df$var05_MiXCR_clonality_group[which(df$var04_MiXCR_clonality >= median_MiXCR)] = "UP"
df$var05_MiXCR_clonality_group[which(df$var04_MiXCR_clonality < median_MiXCR)] = "DOWN"
df$var05_MiXCR_clonality_group = factor(df$var05_MiXCR_clonality_group, levels = c("DOWN", "UP"))
# check
table(df$var05_MiXCR_clonality_group)

# Adaptive TCR clonality
df$var06_Adaptive_TCR_clonality = TCR_Overview$productive_clonality[match(df$Patient_ID, TCR_Overview$Patient_ID)]
df$var07_Adaptive_TCR_clonality_group = NA
median_TCR_Adaptive = median(TCR_Overview$productive_clonality, na.rm = TRUE)
df$var07_Adaptive_TCR_clonality_group[which(df$var06_Adaptive_TCR_clonality >= median_TCR_Adaptive)] = "UP"
df$var07_Adaptive_TCR_clonality_group[which(df$var06_Adaptive_TCR_clonality < median_TCR_Adaptive)] = "DOWN"
df$var07_Adaptive_TCR_clonality_group = factor(df$var07_Adaptive_TCR_clonality_group, levels = c("DOWN", "UP"))
# check
table(df$var07_Adaptive_TCR_clonality_group)

# GIE
df$var08_GIE_value = frequency_df$GIE[match(df$Patient_ID, frequency_df$Patient_ID)]
df$var09_GIE_category = frequency_df$GIE_cat[match(df$Patient_ID, frequency_df$Patient_ID)]
df$var09_GIE_category = factor(df$var09_GIE_category, levels = c("non GIE", "GIE"))
#check
table(df$var09_GIE_category, exclude = NULL)


df$var10_IES = IES_df$IES[match(df$Patient_ID, IES_df$Patient_ID)]
df$var10_IES = factor(df$var10_IES, levels = c("IES1", "IES2", "IES3", "IES4"))
df$var10_IES = as.numeric(df$var10_IES)
#df$var10_IES[which(df$var02_ICR_cluster == "ICR Medium")] = "ICR Medium"
#check
table(df$var10_IES)

# WES data
df$var11_Nonsynonymous_Mutation_count = frequency_df$Non_silent_Mutation_frequency[match(df$Patient_ID, frequency_df$Patient_ID)]
df$var12_Hypermutation_status = frequency_df$Mutation_cat[match(df$Patient_ID, frequency_df$Patient_ID)]
df$var12_Hypermutation_status = factor(df$var12_Hypermutation_status, levels = c("nonhypermutated", "hypermutated"))
df$var13_Neoantigen_count = frequency_df$Neoantigen_count[match(df$Patient_ID, frequency_df$Patient_ID)]
df$var14_Ratio_Neoantigen_count_Nonsynonymous_Mut_count = df$var13_Neoantigen_count / df$var11_Nonsynonymous_Mutation_count
# check
table(df$var12_Hypermutation_status, exclude = NULL)

# Microbiome
Microbiome = rbind(Training, Validation)
df$var15_Microbiome_risk_score = Microbiome$prediction[match(df$Patient_ID, Microbiome$Patient_ID)]
df$var16_Microbiome_risk_group = Microbiome$risk_lab[match(df$Patient_ID, Microbiome$Patient_ID)]
df$var16_Microbiome_risk_group = factor(df$var16_Microbiome_risk_group, 
                                           levels = c("High Risk", "Low Risk"))


# Clinical data
df$var17_Stage_ordinal = clinical_data$ajcc_pathologic_tumor_stage[match(df$Patient_ID,
                                                                         clinical_data$Patient_ID)]
df$var18_Stage_categorical = NA
df$var18_Stage_categorical[which(df$var17_Stage_ordinal %in% c("1", "2"))] = "Stage I-II"
df$var18_Stage_categorical[which(df$var17_Stage_ordinal %in% c("3", "4"))] = "Stage III-IV"
df$var18_Stage_categorical = factor(df$var18_Stage_categorical, levels = c("Stage I-II", "Stage III-IV"))
df$var17_Stage_ordinal = as.numeric(df$var17_Stage_ordinal)

# Add MSI
MANTIS = MANTIS[which(MANTIS$Tissue == "T"),]
MANTIS = MANTIS[which(MANTIS$Patient_ID %in% frequency_df$Patient_ID),]
df$var19_MSI_MANTIS = MANTIS$MSI[match(df$Patient_ID, MANTIS$Patient_ID)]
df$var19_MSI_MANTIS = factor(df$var19_MSI_MANTIS, levels = c("MSS", "MSI-H"))

# ICR ordinal
df$var20_ICR_ordinal = df$var02_ICR_cluster
df$var20_ICR_ordinal = factor(df$var20_ICR_ordinal, levels = c("ICR Low", "ICR Medium", "ICR High"))
df$var20_ICR_ordinal = as.numeric(df$var20_ICR_ordinal)

# More clinical data
df$var21_age = as.numeric(clinical_data$age_at_initial_pathologic_diagnosis[match(df$Patient_ID, clinical_data$Patient_ID)])
df$var22_gender = clinical_data$gender[match(df$Patient_ID, clinical_data$Patient_ID)]
df$var22_gender = factor(df$var22_gender, levels = c("MALE", "FEMALE"))

# CMS
df$var23_CMS = Rfcms$RF.predictedCMS[match(df$Patient_ID, substring(rownames(Rfcms), 1, 3))]
df$var23_CMS[which(is.na(df$var23_CMS))] = "mixed"
df$var23_CMS = factor(df$var23_CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "mixed"))
levels(df$var23_CMS) = c("Other", "Other", "Other", "CMS4", "Other")

dir.create("./Analysis/Multiomics_Survival_Model", showWarnings = FALSE)
dir.create("./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction", showWarnings = FALSE)
save(df, file = "./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.Rdata")
write.csv(df, file = "./Analysis/Multiomics_Survival_Model/001_Input_Data_For_Survival_Prediction/Input_Data_For_Survival_Prediction.csv", row.names = FALSE)

