
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggpubr", "ggplot2"))

# Set parameters


# Load data
load(paste0("./Analysis/WES/020_Mut_Load_and_Neoantigen/Nonsilent_mutation_frequency_and_filtered_Neoantigen_count.Rdata"))
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")

# Hypermutation status
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

# Ratio calculation
frequency_df$Ratio = frequency_df$Neoantigen_count / frequency_df$Non_silent_Mutation_frequency

#### Alternative #####

fit = lm(Neoantigen_count ~ Non_silent_Mutation_frequency, data = frequency_df)
summary(fit)
print(fit)
# Formula: Neoantigen_count = -2.38770 + (0.09171  * Non_silent_Mutation_frequency)

frequency_df$Expected_neoantigen_count = -2.38770 + (0.09171 * frequency_df$Non_silent_Mutation_frequency)
  
frequency_df$GIE = (frequency_df$Neoantigen_count / frequency_df$Expected_neoantigen_count)
frequency_df$GIE_cat = NA
frequency_df$GIE_cat[which(frequency_df$GIE < 1)] = "GIE"
frequency_df$GIE_cat[which(frequency_df$GIE >= 1)] = "non GIE"

plot = ggplot(data = frequency_df, aes(x = Non_silent_Mutation_frequency, y = Neoantigen_count)) +
  stat_cor(method = "pearson", size = 6) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        aspect.ratio = 1/1) +
  #geom_smooth(method="lm") +
  geom_point(aes(color = GIE_cat)) +
  scale_color_manual(values = c("orange", "darkblue"))

dir.create("./Figures/WES/101_AC-ICAM_Calculate_GIE_Clean", showWarnings = FALSE)
png("./Figures/WES/101_AC-ICAM_Calculate_GIE_Clean/101_Scatterplot_Non_silent_Mutation_frequency_by_Neoantigen_count.png",
    res = 600, units = "in", width = 4, height = 4)
plot(plot)
dev.off()

dir.create("./Analysis/WES/101_AC-ICAM_Calculate_GIE_Clean", showWarnings = FALSE)
save(frequency_df, file = "./Analysis/WES/101_AC-ICAM_Calculate_GIE_Clean/frequency_df_GIE_AC_ICAM.Rdata")
