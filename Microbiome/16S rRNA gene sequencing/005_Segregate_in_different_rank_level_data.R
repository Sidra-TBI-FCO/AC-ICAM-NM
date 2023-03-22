
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# Set parameters
Type = "Relative"  # "Normalized" or "Relative"
Mode = "all_samples" # "paired" or "all_tumor" or " "all_samples"

# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type,"_abundancies/clean_", Mode, "_V3.2_sil_OTU_table.Rdata"))

tax_table_dat[which(is.na(tax_table_dat))] = "Unknown"
tax_table_dat = data.frame(tax_table_dat)
table(tax_table_dat[, "Kingdom"], exclude = NULL)
table(tax_table_dat[, "Phylum"], exclude = NULL)
table(tax_table_dat[, "Class"], exclude = NULL) # <NA> 57
table(tax_table_dat[, "Order"], exclude = NULL) # <NA> 129
table(tax_table_dat[, "Family"], exclude = NULL) # <NA> 476
table(tax_table_dat[, "Genus"], exclude = NULL) # <NA> 479
table(tax_table_dat[, "Species"], exclude = NULL) # <NA> 10799

# Class full
tax_table_dat$Class_full = paste(tax_table_dat$Kingdom, tax_table_dat$Phylum,
                                 tax_table_dat$Class)
tax_table_dat$Class_full[which(tax_table_dat$Class == "Unknown")] = "Unknown" 

# Order full
tax_table_dat$Order_full = paste(tax_table_dat$Kingdom, tax_table_dat$Phylum,
                                 tax_table_dat$Class, tax_table_dat$Order)
tax_table_dat$Order_full[which(tax_table_dat$Order == "Unknown")] = "Unknown" 

# Family full
tax_table_dat$Family_full = paste(tax_table_dat$Kingdom, tax_table_dat$Phylum,
                                  tax_table_dat$Class, tax_table_dat$Order,
                                  tax_table_dat$Family)
tax_table_dat$Family_full[which(tax_table_dat$Family == "Unknown")] = "Unknown" 

# Genus full
tax_table_dat$Genus_full = paste(tax_table_dat$Kingdom, tax_table_dat$Phylum,
                                  tax_table_dat$Class, tax_table_dat$Order,
                                  tax_table_dat$Family, tax_table_dat$Genus)
tax_table_dat$Genus_full[which(tax_table_dat$Genus == "Unknown")] = "Unknown"

# Species full
tax_table_dat$Species_full = paste(tax_table_dat$Kingdom, tax_table_dat$Phylum,
                                 tax_table_dat$Class, tax_table_dat$Order,
                                 tax_table_dat$Family, tax_table_dat$Genus,
                                 tax_table_dat$Species)
tax_table_dat$Species_full[which(tax_table_dat$Species == "Unknown")] = "Unknown" 


dat = as.data.frame(OTU_table_dat)

Ranks = colnames(tax_table_dat)
N.ranks = length(Ranks)

i = 8
for (i in 1:N.ranks){
  rank = Ranks[i]
  dat$Rank = tax_table_dat[, rank][match(rownames(dat), rownames(tax_table_dat))]
  rank_abundance = aggregate(.~Rank, dat, FUN=sum)
  rownames(rank_abundance) = rank_abundance$Rank
  rank_abundance$Rank = NULL
  assign(paste0(rank, "_abundance"), rank_abundance)
}

save(Kingdom_abundance, Phylum_abundance, Class_abundance, Order_abundance, Family_full_abundance,
     Genus_full_abundance, Species_full_abundance, file = paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type,"_abundancies/clean_", Mode, "_rank_level_abundancies.Rdata"))

