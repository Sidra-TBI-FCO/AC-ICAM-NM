
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("GenVisR", "dplyr", "ggplot2", "maftools", "ComplexHeatmap",
                      "data.table", "tidyr", "stringr", "cowplot", "tidyverse")
ipak(required.packages)

# Load data
load("./Processed_Data/WES/MAF/MAF_mutect_all_version_3.Rdata")
load("./Processed_Data/WES/MAF/MAF_strelka2_all_version_3.Rdata")

# test = aggregate(data = mutect_MAF_df, t_depth ~ Sample_ID, FUN = sum)

# Filter out SNP and Low complexity variants for strelka2
maf_strelka_INDELS = strelka2_MAF_df %>% filter(Variant_Type != "SNP") 
table(maf_strelka_INDELS$Variant_Classification) # Print table for reporting
nrow(maf_strelka_INDELS) # number of indel variants # 140491

maf_strelka_INDELS = maf_strelka_INDELS %>% filter(!str_detect(DOMAINS, "Low_complexity"))
table(maf_strelka_INDELS$Variant_Classification) # Print table for reporting
nrow(maf_strelka_INDELS) # number of indel variants # 133679
unique(mutect_MAF_df$t_ref_count) # check if all can be converted with as.numeric
mutect_MAF_df$t_ref_count = as.numeric(mutect_MAF_df$t_ref_count)
unique(mutect_MAF_df$n_ref_count) # check if all can be converted with as.numeric
mutect_MAF_df$n_ref_count = as.numeric(mutect_MAF_df$n_ref_count)

unique(maf_strelka_INDELS$t_ref_count) # all values are . , don't convert to numeric!
maf_strelka_INDELS$t_ref_count = as.character(maf_strelka_INDELS$t_ref_count)
unique(maf_strelka_INDELS$n_ref_count) # all values are . , don't convert to numeric!
maf_strelka_INDELS$n_ref_count = as.character(maf_strelka_INDELS$n_ref_count)

# Combine strelka2 results with mutect results by rbind
finalMaf = rbind(maf_strelka_INDELS, mutect_MAF_df)

# VAF filter
unique(finalMaf$t_alt_count) # check if all can be converted with as.numeric
finalMaf$t_alt_count = as.numeric(finalMaf$t_alt_count)
unique(finalMaf$t_depth)
finalMaf$t_vaf = round(finalMaf$t_alt_count / finalMaf$t_depth, 3)

unique(finalMaf$n_alt_count) # check if all can be converted with as.numeric
finalMaf$n_alt_count = as.numeric(finalMaf$n_alt_count)
unique(finalMaf$n_depth)
finalMaf$n_vaf = round(finalMaf$n_alt_count / finalMaf$n_depth, 3)

finalMaf = finalMaf[which(finalMaf$t_vaf > 0.05),] # Min allele fraction > 5%
# This was not applied: finalMaf = finalMaf[-which(finalMaf$n_alt_count >1),] # remove variants with n_alt_count > 1

if(min(finalMaf$t_depth) < 3){
  finalMaf = finalMaf[-which(finalMaf$t_depth < 3), ] # remove variant with t_depth < 3
} 

table(finalMaf$Variant_Classification) # for reporting
nrow(finalMaf)  # for reporting  # 578401

# Filter ExAC
finalMaf = finalMaf  %>%  
  filter (ExAC_AF < 0.01 | is.na(ExAC_AF))
nrow(finalMaf) # for reporting
table(finalMaf$Variant_Classification) # for reporting

# Filter MAF file for non-synonymous variants
finalMafFiltered = finalMaf  %>%
  filter (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", 
                                        "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))
nrow(finalMafFiltered) # for reporting
table(finalMafFiltered$Variant_Classification) # for reporting

save(finalMafFiltered, file = "./Processed_Data/WES/MAF/finalMafFiltered_nonsynonymous_filter_all_samples.Rdata")

# Filter not tolerated or not benign
finalMafFiltered = finalMafFiltered  %>% 
  filter (!str_detect(SIFT , "tolerated") | !str_detect(PolyPhen, "benign"))
nrow(finalMafFiltered) # for reporting

# Just for reporting
genes = finalMafFiltered[, c("Hugo_Symbol"), drop = FALSE]
artifacts = c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG",paste0("OR",1:9))
genes = genes[-which(substring(genes$Hugo_Symbol, 1, 3) %in% artifacts), , drop = FALSE]
artifacts2 = c("PLIN","CELA","SRA1")
genes = genes[-which(substring(genes$Hugo_Symbol, 1, 4) %in% artifacts2), , drop = FALSE]
artifacts_genes = c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")
genes = genes[-which(genes$Hugo_Symbol %in% artifacts_genes),]
length(genes)

save(finalMaf, file = "./Processed_Data/WES/MAF/finalMaf_all_samples.Rdata")

save(finalMafFiltered, file = "./Processed_Data/WES/MAF/finalMafFiltered_all_samples.Rdata")

