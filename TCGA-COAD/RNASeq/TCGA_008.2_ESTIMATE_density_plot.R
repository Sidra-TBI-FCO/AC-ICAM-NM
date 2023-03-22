
# Setup environment
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("corrplot", "stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
variable_of_interest = "ESTIMATEScore" #"ESTIMATEScore" #"ImmuneScore" #"ESTIMATEScore"
intersect = ""

# Create directories
dir.create("./Figures/Trimmed_p/ESTIMATE", showWarnings = FALSE)

# Load data
load("./Analysis/008_ESTIMATE/Biolinks_TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")
ESTIMATE_Bio = ESTIMATE
load("../NGS_Data/Analysis/Trimmed_p/ESTIMATE/TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")
ESTIMATE_Assembler = ESTIMATE_tcga

ESTIMATE_Bio = as.data.frame(ESTIMATE_Bio)
ESTIMATE_Assembler = as.data.frame(ESTIMATE_Assembler)
rownames(ESTIMATE_Assembler) = gsub("\\.", "-", rownames(ESTIMATE_Assembler))

ESTIMATE_Bio = ESTIMATE_Bio[-which(rownames(ESTIMATE_Bio) %in% rownames(ESTIMATE_Assembler)),]

ESTIMATE_Bio$cohort = "Biolinks"
ESTIMATE_Assembler$cohort = "TCGA_Assembler"
df_plot = rbind(ESTIMATE_Bio, ESTIMATE_Assembler)

plot = ggplot(df_plot, aes(x = df_plot[,variable_of_interest], color = cohort)) +
  geom_density(aes(fill = cohort), alpha = .4) +
  xlab(paste0(variable_of_interest)) +
  scale_color_manual(values = c("#CC99FF", "lightblue")) +
  scale_fill_manual(values = c("#CC99FF", "lightblue")) +
  theme_bw() +
  xlab("ESTIMATE score") +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15))

dir.create("./Figures/ESTIMATE", showWarnings = FALSE)
png(paste0("./Figures/ESTIMATE/", "Density_plot_", variable_of_interest, "_Biolinks_vs_Assembler.png"), res = 600, height = 4, width = 5, units = "in")
plot(plot)
dev.off()
