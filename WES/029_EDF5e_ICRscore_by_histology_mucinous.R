
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

required.packages = c("ggplot2", "ggpubr")
ipak(required.packages)

# Load data
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")

BRCA2_mutated_samples = finalMafFiltered$Patient_ID[which(finalMafFiltered$Hugo_Symbol == "BRCA2")]
HR_mutated_samples = finalMafFiltered$Patient_ID[which(finalMafFiltered$Hugo_Symbol %in% c("BRCA2", "FANCA", "BRCA1"))]

# Prepare data
df = frequency_df
df$Mutation_cat = NA
df$Mutation_cat[which(df$Nonsilent_mutational_burden_per_Mb <= 12)] = "nonhypermutated"
df$Mutation_cat[which(df$Nonsilent_mutational_burden_per_Mb > 12)] = "hypermutated"
table(df$Mutation_cat)

df$ICRscore = table_cluster_assignment$ICRscore[match(df$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
df$Histology = clinical_data$Tumor_morphology[match(df$Patient_ID, clinical_data$Patient_ID)]
df$Histology_simple = df$Histology
df$Histology_simple[which(df$Histology == "mucineus adenocarcinoom")] = "Mucinous"
df$Histology_simple[-which(df$Histology == "mucineus adenocarcinoom")] = "Others"

df$BRCA2_status = "WT"
df$BRCA2_status[which(df$Patient_ID %in% BRCA2_mutated_samples)] = "BRCA2_MUT"

df$HR_status = "WT"
df$HR_status[which(df$Patient_ID %in% HR_mutated_samples)] = "HR_MUT"

df = df[which(df$Mutation_cat == "hypermutated"),]


# Plot
plot = ggplot(df, aes(x = HR_status, y = ICRscore, fill = Histology_simple)) +
  scale_fill_manual(values = c("Mucinous" = "darkorange", "Others" = "darkgrey")) +
  facet_grid(~Histology_simple) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, aes(color = HR_status)) +
  xlab("") +
  ylab("") +
  #stat_compare_means(comparisons = list(c("HR_MUT", "WT")),
   #                  label = "p.signif", label.y =  10) +
  theme_bw() +
  ylim(4, 11) +
  ylab("ICR score") +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black"))

dir.create("./Figures/WES/029_ICRscore_by_mucinous_subtype", showWarnings = FALSE)  
pdf("./Figures/WES/029_ICRscore_by_mucinous_subtype/EDF5e_Boxplot_ICRscore_by_Mucinous_and_others_in_hypermutated.pdf",
    #width = 2.4, height = 3, pointsize = 12)
    width = 3.5, height = 3)
plot(plot)
dev.off()

chisq.test(table(df$Histology_simple, df$HR_status))
table(df$Histology_simple, df$HR_status)

chisq.test(table(df$Histology_simple, df$BRCA2_status))
table(df$Histology_simple, df$BRCA2_status)
