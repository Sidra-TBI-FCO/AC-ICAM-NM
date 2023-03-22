
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("GenVisR", "dplyr", "ggplot2", "maftools", "ComplexHeatmap",
                      "data.table", "tidyr", "stringr", "cowplot", "tidyverse")
ipak(required.packages)

# Download package TCGAmutations
devtools::install_github(repo = "PoisonAlien/TCGAmutations")
library("TCGAmutations")

# Test tcgaCompare function with tcga COAD from MC3 file
tcga_available()
tcga_load(study = "COAD", source = "MC3")
TCGA_COAD_MC3_df = tcga_coad_mc3@data
table(TCGA_COAD_MC3_df$Variant_Classification)
table(TCGA_COAD_MC3_df$ExAC_AF)
variants_per_sample = tcga_coad_mc3@variants.per.sample
max(variants_per_sample$Variants)
min(variants_per_sample$Variants)
coad = read.maf(maf = TCGA_COAD_MC3_df)

dir.create("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/WES", showWarnings = FALSE)
dir.create("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/WES/MAF", showWarnings = FALSE)
save(TCGA_COAD_MC3_df, file = "../NGS_Data_TCGA_COAD_Jessica/Processed_Data/WES/MAF/TCGAmutation_MAF_TCGA_COAD.Rdata")

load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/WES/MAF/TCGAmutation_MAF_TCGA_COAD.Rdata")

png("./Figures/WES/002_tcgaCompare_TMB/Test-TCGA-COAD-tcgaCompare_TMB.png", res = 600, 
    height = 6, width = 12, units = "in")
coad.mutload = tcgaCompare(maf = coad, cohortName = 'COAD-test')
dev.off()

png("./Figures/WES/002_tcgaCompare_TMB/Test-TCGA-COAD-tcgaCompare_TMB_per_Mb_capture_size_50.png", res = 600, 
    height = 6, width = 12, units = "in")
coad.mutload = tcgaCompare(maf = coad, cohortName = 'COAD-test', tcga_capture_size = 50, capture_size = 50)
dev.off()

png("./Figures/WES/002_tcgaCompare_TMB/Test-TCGA-COAD-tcgaCompare_TMB_per_Mb_capture_size_40.png", res = 600, 
    height = 6, width = 12, units = "in")
coad.mutload = tcgaCompare(maf = coad, cohortName = 'COAD-test', tcga_capture_size = 40, capture_size = 40)
dev.off()

# Load data
load("./Processed_Data/WES/MAF/finalMaf_clean_primary_tumor_samples.Rdata")
finalMaf$Tumor_Sample_Barcode = as.character(finalMaf$Tumor_Sample_Barcode)

finalMaf = finalMaf %>%
  filter (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                        "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
                                        "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))

table(finalMaf$Variant_Classification)

silu= read.maf(maf = finalMaf)

dir.create("./Figures/WES/002_tcgaCompare_TMB", showWarnings = FALSE)

pdf("./Figures/WES/002_tcgaCompare_TMB/v6_AC_ICAM_tcgaCompare_TMB_per_Mb_capture_size_40.pdf", 
    height = 4, width = 8)
silu.mutload = tcgaCompare(maf = silu, cohortName = 'AC-ICAM', tcga_capture_size = 40, capture_size = 40)
dev.off()
