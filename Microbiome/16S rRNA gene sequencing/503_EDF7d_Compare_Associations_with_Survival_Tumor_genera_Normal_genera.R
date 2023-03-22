
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak("VennDiagram")

# Set parameters
Survival = "PFS"

# Load data
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
tumor_genera = rownames(Genus_full_abundance)
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_N_246_samples_based_on_normal.Rdata")
normal_genera = rownames(Genus_full_abundance)

load("./Analysis/Microbiome/021_KM_by_median/PFS_KM_relative_abundance_by_median.Rdata")
PFS_T = results
load("./Analysis/Microbiome/021_KM_by_median/OS_KM_relative_abundance_by_median.Rdata")
OS_T = results

load("./Analysis/Microbiome/021_KM_by_median/PFS_KM_relative_abundance_by_median_N.Rdata")
PFS_N = results
load("./Analysis/Microbiome/021_KM_by_median/OS_KM_relative_abundance_by_median_N.Rdata")
OS_N = results

# Analysis
analysis_df = data.frame(Tissue = c("T", "T", "N", "N"), Survival = c("OS", "PFS", "OS", "PFS"),
                         Number_significant = NA, Number_FDR = NA)

# Tumor
analysis_df[which(analysis_df$Tissue == "T" & analysis_df$Survival == "OS"), "Number_significant"] = length(which(OS_T$p_val < 0.05))
analysis_df[which(analysis_df$Tissue == "T" & analysis_df$Survival == "PFS"), "Number_significant"] = length(which(PFS_T$p_val < 0.05))
analysis_df[which(analysis_df$Tissue == "T" & analysis_df$Survival == "OS"), "Number_FDR"] = length(which(OS_T$FDR < 0.1))
analysis_df[which(analysis_df$Tissue == "T" & analysis_df$Survival == "PFS"), "Number_FDR"] = length(which(PFS_T$FDR < 0.1))

# Normal
analysis_df[which(analysis_df$Tissue == "N" & analysis_df$Survival == "OS"), "Number_significant"] = length(which(OS_N$p_val < 0.05))
analysis_df[which(analysis_df$Tissue == "N" & analysis_df$Survival == "PFS"), "Number_significant"] = length(which(PFS_N$p_val < 0.05))
analysis_df[which(analysis_df$Tissue == "N" & analysis_df$Survival == "OS"), "Number_FDR"] = length(which(OS_N$FDR < 0.05))
analysis_df[which(analysis_df$Tissue == "N" & analysis_df$Survival == "PFS"), "Number_FDR"] = length(which(PFS_N$FDR < 0.05))

class(analysis_df$Number_significant) = "numeric"

OS_T_signif = OS_T$Name[which(OS_T$p_val < 0.05)]
OS_N_signif = OS_N$Name[which(OS_N$p_val < 0.05)]

PFS_T_signif = PFS_T$Name[which(PFS_T$p_val < 0.05)]
PFS_N_signif = PFS_N$Name[which(PFS_N$p_val < 0.05)]

library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(x = list(OS_N_signif, OS_T_signif, PFS_T_signif, PFS_N_signif),
             category.names = c("OS_N_signif", "OS_T_signif", "PFS_T_signif", "PFS_N_signif"),
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             filename = "./Figures/Microbiome/New_Compare_Association_with_Survival_By_Median/VennDiagram_Genera.png",
             output = TRUE,
             # Set names
             cat.fontfamily = "sans")


#analysis_df = analysis_df[which(analysis_df$Survival == Survival),]

plot = ggplot(analysis_df, aes(x = Tissue, y = Number_significant, fill = Tissue)) +
  scale_fill_manual(values = c("T" = "#27A4DE", "N" = "#CC578F")) +
  ylab("") +
  facet_grid(.~Survival) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size = 12))

dir.create("./Figures/Microbiome/New_Compare_Association_with_Survival_By_Median", showWarnings = FALSE)
png(paste0("./Figures/Microbiome/New_Compare_Association_with_Survival_By_Median/Compare_Association_with_Survival_", Survival,"By_Median.png"),
    res = 600, units = "in", width = 3, height = 3)
plot(plot)
dev.off()

pdf(paste0("./Figures/Microbiome/New_Compare_Association_with_Survival_By_Median/Compare_Association_with_Survival_Both_PFS_and_OS_", Survival,"By_Median.pdf"),
    width = 4, height = 3)
plot(plot)
dev.off()

OS_T = OS_T[which(OS_T$Name %in% tumor_genera),]
PFS_T = PFS_T[which(PFS_T$Name %in% tumor_genera),]

OS_N = OS_N[which(OS_N$Name %in% normal_genera),]
PFS_N = PFS_N[which(PFS_N$Name %in% normal_genera),]

# Calculate FDR
OS_T$FDR = p.adjust(OS_T$p_val, method = "BH", n = nrow(OS_T))
OS_T = OS_T[, c(1, 3, 4, 5, 2, 6)]
OS_N$FDR = p.adjust(OS_N$p_val, method = "BH", n = nrow(OS_N))
OS_N = OS_N[, c(1, 3, 4, 5, 2, 6)]
PFS_T$FDR = p.adjust(PFS_T$p_val, method = "BH", n = nrow(PFS_T))
PFS_T = PFS_T[, c(1, 3, 4, 5, 2, 6)]
PFS_N$FDR = p.adjust(PFS_N$p_val, method = "BH", n = nrow(PFS_N))
PFS_N = PFS_N[, c(1, 3, 4, 5, 2, 6)]

dir.create("./Analysis/Microbiome/New_Associations_with_Survival", showWarnings = FALSE)
write.csv(OS_T, file = "./Analysis/Microbiome/New_Associations_with_Survival/Coxph_OS_based_on_median_T.csv",
          row.names = FALSE)
write.csv(OS_N, file = "./Analysis/Microbiome/New_Associations_with_Survival/Coxph_OS_based_on_median_N.csv",
          row.names = FALSE)
write.csv(PFS_T, file = "./Analysis/Microbiome/New_Associations_with_Survival/Coxph_PFS_based_on_median_T.csv",
          row.names = FALSE)
write.csv(PFS_N, file = "./Analysis/Microbiome/New_Associations_with_Survival/Coxph_PFS_based_on_median_N.csv",
          row.names = FALSE)
