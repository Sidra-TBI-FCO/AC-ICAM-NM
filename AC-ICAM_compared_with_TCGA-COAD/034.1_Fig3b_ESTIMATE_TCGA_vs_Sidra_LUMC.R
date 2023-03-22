
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "ESTIMATEScore" #"ImmuneScore" #"StromalScore" #"ESTIMATEScore"

# Create directories
dir.create("./Figures/Trimmed_p/034.1_ESTIMATE_TCGA_vs_Sidra_LUMC", showWarnings = FALSE)

# Load data
load(paste0("./Analysis/Trimmed_p/034_ESTIMATE_Merged_TCGA_SILU/original_ESTIMATE.score.Rdata"))

ESTIMATE = as.data.frame(ESTIMATE)

ESTIMATE$cohort = "AC-ICAM"
ESTIMATE$cohort[grepl("TCGA", rownames(ESTIMATE))] = "TCGA-COAD"
df_plot = ESTIMATE

my_comparisons = list(c("SILU-COAD", "TCGA-COAD"))
plot = ggplot(df_plot, aes(x = cohort, y = ESTIMATEScore)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
  geom_jitter(size = 0.2, width = 0.2) +
  theme_bw() +
  ylim(min(df_plot$ESTIMATEScore)-100, 5000) +
  ylab("ESTIMATE score") +
  scale_fill_manual(values = c("#F2A481", "#BCEDED")) +
  xlab("") +
  theme(axis.text.x = element_text(color = "black", size = 15,
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        legend.position = "none") 

#plot = ggboxplot(df_plot, x = "cohort", y = variable_of_interest,
 #                color = "cohort", palette = c("Sidra-LUMC" = "orange", "TCGA" = "#006DFF"),
   #              add = "jitter", outlier.shape = NA) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test")

png(paste0("./Figures/Trimmed_p/034.1_ESTIMATE_TCGA_vs_Sidra_LUMC/v4_", variable_of_interest, "_JSREP_vs_TCGA.png"), 
    res = 600, height = 4, width = 3, units = "in")
plot(plot)
dev.off()

pdf(paste0("./Figures/Trimmed_p/034.1_ESTIMATE_TCGA_vs_Sidra_LUMC/v4_", variable_of_interest, "_AC_ICAM_vs_TCGA.pdf"), 
    height = 4, width = 3)
plot(plot)
dev.off()

test = t.test(df_plot$ESTIMATEScore[which(df_plot$cohort == "AC-ICAM")], df_plot$ESTIMATEScore[which(df_plot$cohort == "TCGA-COAD")],
              alternative = "two.sided", paired = FALSE)
test$p.value


