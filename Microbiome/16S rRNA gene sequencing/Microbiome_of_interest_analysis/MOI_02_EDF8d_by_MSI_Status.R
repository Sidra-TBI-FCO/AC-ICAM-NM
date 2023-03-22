

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2"))

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")

# Set parameters 
# Genus should be one of the rownames of Genus_full_abundance
# rownames(Genus_full_abundance)[grep("Fuso", rownames(Genus_full_abundance))
Genus = "D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"
  #"D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 1"
#"D_0__Bacteria D_1__Fusobacteria D_2__Fusobacteriia D_3__Fusobacteriales D_4__Fusobacteriaceae D_5__Fusobacterium"

dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots", showWarnings = FALSE)
dir.create(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_02_Genus_by_MSI"), showWarnings = FALSE)

# Analysis
frequency_df$Mutation_cat = NA
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
frequency_df$Mutation_cat[which(frequency_df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"

frequency_df = frequency_df[which(frequency_df$Patient_ID %in% 
                                    substring(colnames(Genus_full_abundance), 1, 3)),]

df_plot = frequency_df

df_plot$MSI_Status = MANTIS$MSI[match(df_plot$Patient_ID, MANTIS$Patient_ID)]

df_plot$Abundance = Genus_full_abundance[Genus,][match(df_plot$Patient_ID,
                                                       substring(colnames(Genus_full_abundance), 1, 3))]

df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4

df_plot$Mutation_cat = factor(df_plot$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))
table(df_plot$Mutation_cat)

df_plot$MSI_Status = factor(df_plot$MSI_Status, levels = c("MSS", "MSI-H"))

plot = ggplot(df_plot, aes(x = MSI_Status, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, aes(fill=MSI_Status)) +
  scale_fill_manual(values = c("MSS" = "yellow", "MSI-H" = "purple")) +
  geom_jitter(width = 0.2, size = 0.8) +
  #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
  theme_bw() +
  scale_y_log10() +
  ylab("") +
  #ylab(paste0(gsub(".*\\D5__", "", Genus), "\n","Relative abundance")) +
  xlab("") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none")

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_02_Genus_by_MSI/EDF8d_1_", Genus, 
           "_by_MSI_Status.pdf"), width = 3, height = 4)
plot(plot)
dev.off()

plot = ggplot(df_plot, aes(x = Mutation_cat, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, aes(fill=Mutation_cat)) +
  scale_fill_manual(values = c("nonhypermutated" = "#A7EABD", "hypermutated" = "#EAAED0")) +
  geom_jitter(width = 0.2, size = 0.8) +
  #scale_color_manual(values = c("nonhypermutated" = "darkgreen", "hypermutated" = "orchid")) +
  theme_bw() +
  scale_y_log10() +
  ylab("") +
  #ylab(paste0(gsub(".*\\D5__", "", Genus), "\n","Relative abundance")) +
  xlab("") +
  theme(axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.position = "none")

pdf(paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_02_Genus_by_MSI/EDF8d_2_", Genus, 
           "_by_Hypermutation_Status.pdf"), width = 3, height = 4)
plot(plot)
dev.off()


# run stats without running line 41  df_plot$Abundance[which(df_plot$Abundance == 0)] = 1e-4
t.test(df_plot$Abundance[which(df_plot$MSI_Status == "MSI-H")],
       df_plot$Abundance[which(df_plot$MSI_Status == "MSS")])


t.test(df_plot$Abundance[which(df_plot$Mutation_cat == "hypermutated")],
       df_plot$Abundance[which(df_plot$Mutation_cat == "nonhypermutated")])
