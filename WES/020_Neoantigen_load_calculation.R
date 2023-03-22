
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("ggplot2", "stringr", "ggpubr")
ipak(required.packages)

# Load data
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
neoantigens_full = read.table(file = "./Processed_Data/WES/Neoantigens_pvactools/Filtered.tsv", sep = "\t", stringsAsFactors = FALSE,
                              header = TRUE)
neoantigens_500 = read.table(file = "./Processed_Data/WES/Neoantigens_pvactools/Filtered_500.tsv", sep = "\t", stringsAsFactors = FALSE,
                             header = TRUE)


# Only mutant peptides with binding affinity IC50 < 500 nM for the restricted HLA-I subtype 
# and binding affinity of the corresponding wild-type peptide IC50 > 500 nM for all HLA-I subtypes were retained as neoantigens
# Ref: https://www.nature.com/articles/s42003-019-0369-7

# Check filtering
neoantigens_NA = neoantigens_full[which(is.na(neoantigens_full$Median.WT.Score)),]
neoantigens_subset = neoantigens_full[which(neoantigens_full$Median.WT.Score > 500), ]
#nrow(neoantigens_subset) == nrow(neoantigens_500) #correct filtering!

min(neoantigens_500$Median.MT.Score) # 2.16
max(neoantigens_500$Median.MT.Score) # 499.972 # Confirmed

# Count number of neoantigens per sample
# Unfiltered
neoantigen_df = data.frame(table(neoantigens_full$Sample.Name))
colnames(neoantigen_df) = c("Sample_ID", "Neoantigen_count_unfiltered")
neoantigen_df$Sample_ID = str_pad(neoantigen_df$Sample_ID, pad = "0", 4)
neoantigen_df$Patient_ID = substring(neoantigen_df$Sample_ID, 1, 3)
neoantigen_df$Tissue_type = substring(neoantigen_df$Sample_ID, 4, nchar(neoantigen_df$Sample_ID))
neoantigen_df = neoantigen_df[, c(1, 3, 4, 2)]

# Filtered
neoantigen_df_filtered = data.frame(table(neoantigens_500$Sample.Name))
colnames(neoantigen_df_filtered) = c("Sample_ID", "Neoantigen_count")
neoantigen_df_filtered$Sample_ID = str_pad(neoantigen_df_filtered$Sample_ID, pad = "0", 4)

# Find out which sample has a filtered neoantigen load of zero
neoantigen_df$Sample_ID[-which(neoantigen_df$Sample_ID %in% neoantigen_df_filtered$Sample_ID)]
# Check in unfiltered file
check = neoantigens_full[which(neoantigens_full$Sample.Name == "277T"),]
max(check$Median.WT.Score) # 492.477
nrow(check)

# Add to neoantigen df
neoantigen_df$Neoantigen_count = neoantigen_df_filtered$Neoantigen_count[match(neoantigen_df$Sample_ID, neoantigen_df_filtered$Sample_ID)]
neoantigen_df$Neoantigen_count[which(neoantigen_df$Sample_ID == "277T")] = 0

# Quick checks

plot = ggplot(data = neoantigen_df, aes(x = Neoantigen_count_unfiltered, y = Neoantigen_count)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Neoantigen count unfiltered") +
  ylab("Neoantigen count \n Median WT score > 500 nM") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

dir.create("./Figures/WES/020_Neoantigen_load_calculation", showWarnings = FALSE)
png(paste0("./Figures/WES/020_Neoantigen_load_calculation/Filtered_vs_Unfiltered_counts.png"),
    res = 600, width = 5, height = 5, units = "in")
plot(plot)
dev.off()

save(neoantigen_df, file = paste0("./Analysis/WES/020_Mut_Load_and_Neoantigen/All_samples_filtered_Neoantigen_count.Rdata"))

neoantigen_df = neoantigen_df[which(neoantigen_df$Tissue_type == "T"),]

frequency_df$Neoantigen_count = neoantigen_df$Neoantigen_count[match(frequency_df$Patient_ID,
                                                                     neoantigen_df$Patient_ID)]



dir.create("./Analysis/WES/020_Mut_Load_and_Neoantigen", showWarnings = FALSE)
save(frequency_df, file = paste0("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata"))


# Scatterplot Mutational load and neoantigen count

plot = ggplot(data = frequency_df, aes(x = Non_silent_Mutation_frequency, y = Neoantigen_count)) +
  geom_point() +
  stat_cor(method = "pearson", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  xlab("Nonsynonymous mutation count") +
  ylab("Neoantigen count") +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

png(paste0("./Figures/WES/020_Neoantigen_load_calculation/Mut_load_versus_Neoantigen_count.png"),
    res = 600, width = 5.2, height = 5, units = "in")
plot(plot)
dev.off()

