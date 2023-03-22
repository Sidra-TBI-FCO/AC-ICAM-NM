
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Rank = "Genus_full" # For "Kingdom" "Phylum" "Class" "Order" "Family_full" "Genus_full" "Species_full" 
                # for Family, Genus, and Species use the full taxa name (some species have the same species name but are in different phylum/class/order)
Tissue = "N"

# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type, "_abundancies/clean_paired_rank_level_abundancies.Rdata"))

rank_abundancies = get(paste0(Rank, "_abundance"))

abundance_T = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "T")]
colnames(abundance_T) = substring(colnames(abundance_T), 1, 3)
abundance_T = t(abundance_T)
abundance_N = rank_abundancies[,which(substring(colnames(rank_abundancies), 4, 4) == "N")]
colnames(abundance_N) = substring(colnames(abundance_N), 1, 3)
abundance_N = t(abundance_N)

abundance = get(paste("abundance_", Tissue, sep = ""))

df = data.frame(Patient_ID = unique(substring(colnames(rank_abundancies), 1, 3)),
                Abundance = NA,
                Category = NA,
                stringsAsFactors = FALSE)
abundance_cat = abundance
abundance_cat[1:nrow(abundance_cat),] = NA

for (i in 1:nrow(rank_abundancies)){
  micr = rownames(rank_abundancies)[i]
  df$Abundance = abundance[, micr][match(df$Patient_ID, rownames(abundance))]
  median = median(df$Abundance)
  if(is.na(median)){next}
  df$Category[which(df$Abundance > median)] = "High"
  df$Category[which(df$Abundance <= median)] = "Low"
  table(df$Category)
  abundance_cat[, micr] = df$Category
}

dir.create("./Analysis/Microbiome/020_Categories_by_median", showWarnings = FALSE)
save(abundance_cat, file= paste0("./Analysis/Microbiome/020_Categories_by_median/020_Categorized_", Type, "_abundance_by_median_", Tissue, "_", Rank,".Rdata"))
