
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak("tidyr")

# load data
load("./Processed_Data/Microbiome/001_data_preparation/datastruct_V3.2_sil_large_df.Rdata")

dim(datastruct)
# 9412980      44


colnames(datastruct)
#[1] "OTU"                                 "Sample"                              "Abundance"                           "BarcodeSequence"                    
#[5] "BarcodeSequencei5"                   "BarcodeSequencei7"                   "Description"                         "ICR_cluster"                        
#[9] "ICRscore"                            "LinkerPrimerSequence"                "MiSeqRun"                            "Paired"                             
#[13] "Patient_ID"                          "Plate"                               "Plate.number"                        "Primer1"                            
#[17] "Primer2"                             "PrimerSet"                           "RCLinkerPrimerSequence"              "RCReversePrimer"                    
#[21] "Replicate"                           "Sample.ID"                           "SampleCode"                          "SampleName"                         
#[25] "SamplesDescription"                  "SampsID"                             "sample_Species"                      "Subject"                            
#[29] "Tissue"                              "Well"                                "age_at_initial_pathologic_diagnosis" "age_category"                       
#[33] "ajcc_pathologic_tumor_stage"         "gender"                              "histology"                           "primary_tumor_anatomic_site"        
#[37] "primary_tumor_side"                  "Kingdom"                             "Phylum"                              "Class"                              
#[41] "Order"                               "Family"                              "Genus"                               "Species"



# There are 49542 rows without Sample.ID. Don't use Sample.ID as sample identifier (is an old version)
table(datastruct$Sample.ID)
#        002N   002T   004T   005T   007N   009N   009T   010N   010T   013T   016T 017LM1   018N   019T   020N   020T   021T   022N   022T   023N   023T   024T   025N   025T ..etc
# 49542  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514  16514 ...etc


datastruct_no_sample = datastruct[which(datastruct$Sample.ID == ""),]
table(datastruct_no_sample$SampleCode)
## 27       450LM1      450LM2 
## 16514    16514       16514 

datastruct_no_tissue = datastruct[which(datastruct$Tissue == ""),]
table(datastruct_no_tissue$SampleCode)
table(datastruct_no_tissue$SamplesDescription)
## 27       450LM1      450LM2 
## 16514    16514       16514 

################################################################################################
# Use SampleCode instead of Sample.ID
# SampleCode is Patient_ID (is complete, no empty fields)
table(datastruct$SampleCode, exclude = NULL)
table(datastruct$SamplesDescription, exclude = NULL)

# SamplesDescription is Tissue Type (is complete, no empty fields)
#   LM1     LM2       N       T 
# 231196   33028 4359696 4789060 

# Combined
datastruct$SampleCode_new = paste(datastruct$SampleCode, datastruct$SamplesDescription, sep = "")

df_n_rows_per_sample = data.frame(table(datastruct$SampleCode_new))
# 16514 rows per sample, confirmed
# total number of samples = 570

table(datastruct$sample_Species)
# Human 
# 9412980 

table(datastruct$Species, exclude = NULL)

datastruct_small = datastruct[, c("SampleCode_new", "OTU", "Abundance")]
#datastruct_small = datastruct_small[-which(is.na(datastruct_small$Species)),]

df.wide = pivot_wider(datastruct_small, names_from = OTU, values_from = Abundance)
df.wide = data.frame(df.wide)
rownames(df.wide) = df.wide$SampleCode_new
df.wide$SampleCode_new = NULL
OTU_matrix = as.matrix(t(df.wide))
# Values in `Abundance` are not uniquely identified