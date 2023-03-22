
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("reshape2", "ggplot2", "ggpubr"))

# Load data
load("./Analysis/ICR Consensus Clustering/TCGA_GDC_Biolinks_COAD_ICR_cluster_assignment_k2-6.Rdata")
TCGA_COAD = table_cluster_assignment
load("../NGS_data/Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("../NGS_data/Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load(paste0(toolbox.path, "/ICR genes/ICR_genes.RData"))

# Prepare data
dim(RNASeq.QN.LOG2)

table_cluster_assignment$Patient_ID = substring(rownames(table_cluster_assignment), 1, 3)
table_cluster_assignment = table_cluster_assignment[which(table_cluster_assignment$Patient_ID %in% 
                                                            substring(colnames(RNASeq.QN.LOG2), 1, 3)),]

AC_ICAM = table_cluster_assignment

df = data.frame(prop.table(table(AC_ICAM$ICR_HML)))
colnames(df) = c("ICR", "AC_ICAM")
df$TCGA_COAD = NA

df_TCGA = data.frame(prop.table(table(TCGA_COAD$ICR_HML)))
df$TCGA_COAD = df_TCGA$Freq[match(df$ICR, df_TCGA$Var1)]

df_melt = melt(df, id.vars = "ICR")
df_melt$ICR = factor(df_melt$ICR, levels = c("ICR High", "ICR Medium", "ICR Low"))

pdf("../NGS_data/Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/ICR_cluster_assignment_by_cohort.pdf",
    width = 3, height = 3)
ggplot(df_melt, aes(x = variable, y = value, fill = ICR)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ICR High" = "red",
                               "ICR Medium" = "green",
                               "ICR Low" = "blue")) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12)) +
  ylab("Proportion") +
  xlab("")
dev.off()

colnames(AC_ICAM)
colnames(TCGA_COAD)
TCGA_COAD$Patient_ID = substring(rownames(TCGA_COAD), 1, 12)

table_cluster_assignment = rbind(AC_ICAM, TCGA_COAD)

save(table_cluster_assignment, file = "../NGS_data/Analysis/030_Merged_TCGA_COAD_Sidra_LUMC/ICR_clusters_rbind.Rdata")

load("../NGS_data/Processed_Data/RNASeq/030_Merged_TCGA_COAD_Sidra_LUMC/TCGA_COAD_Sidra_LUMC_Filtered_Primary_tumor_EDASEQ_QN_LOG2.Rdata")
colnames(RNASeq.QN.LOG2)
ICR_subset = RNASeq.QN.LOG2[ICR_genes,]

df = data.frame(Sample = colnames(ICR_subset), Cohort = "AC-ICAM", ICR_cluster = NA, ICRscore = colMeans(ICR_subset))

df$Cohort[grep("TCGA", df$Sample)] = "TCGA-COAD"
df$ICR_cluster = table_cluster_assignment$ICR_HML[match(df$Sample, rownames(table_cluster_assignment))]

df$ICR_cluster = factor(df$ICR_cluster, levels = c("ICR Low", "ICR Medium", "ICR High"))
my_comparisons = list(c("TCGA-COAD", "AC-ICAM"))

#df = df[-which(df$ICR_cluster %in% c("ICR High", "ICR Medium")),]
plot = ggplot(df, aes(x = Cohort, y = ICRscore, fill = ICR_cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.8, width = 0.2) +
  scale_fill_manual(values = c("blue", "green", "red")) +
  facet_grid(.~ICR_cluster) +
  xlab("") +
  ylim(3.8, 12) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  theme_bw() +
  theme(axis.title.x = element_text(colour = "black", size = 14),
        axis.title.y = element_text(colour = "black", size = 14),
        axis.text.x = element_text(colour = "black", size = 14,
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 14),
        strip.text = element_text(colour = "black", size = 14),
        strip.background = element_blank(),
        legend.position = "none")

png("../NGS_data/Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/ICRscore_by_cohort.png",
    res = 600, units = "in", width = 6, height = 4)
plot(plot)  
dev.off()

86/348
186/(186+146+107)

levels(df$ICR_cluster) =c("ICR Low", "Other", "Other")

chisq.test(table(df$Cohort, df$ICR_cluster))
df_tcga = df[which(df$Cohort == "TCGA-COAD"),]
df_AC_ICAM = df[which(df$Cohort == "AC-ICAM"),]

cv_tcga <- sd(df_tcga$ICRscore) / mean(df_tcga$ICRscore) * 100

cv_tcga

cv <- sd(df_AC_ICAM$ICRscore) / mean(df_AC_ICAM$ICRscore) * 100
cv



