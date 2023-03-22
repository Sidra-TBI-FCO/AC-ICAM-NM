
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
Source = "ConsensusTME_COAD"

# Load data
results_ICR = read.csv(paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_ICR.csv"), stringsAsFactors = FALSE)
results_templates = read.csv(paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_productive_templates.csv"), stringsAsFactors = FALSE)
#write.csv(results_productive_templates_cor, file = paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_productive_templates_cor.csv"), row.names = FALSE)
results_clonality = read.csv(paste0("./Analysis/TCR/006_T_cell_infiltration_cor/", Source, "results_clonality.csv"), stringsAsFactors = FALSE)

results_templates$Variable = gsub("\\_", " ", results_templates$Variable)
results_clonality$Variable = gsub("\\_", " ", results_clonality$Variable)
cells = results_clonality$Variable[order(results_clonality$cor_coef)]
results_templates$group = "productive_templates"
results_clonality$group = "productive_clonality"
df_plot = rbind(results_templates, results_clonality)
df_plot$Variable = factor(df_plot$Variable, levels = cells)

colors = rep("black", 20)
#colors[grepl("T ", levels(df_plot$Variable))] = "darkblue"

# Plotting
plot = ggplot(df_plot, aes(x= Variable, y = cor_coef, group = group, color = group)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_line() +
  geom_point() +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("darkorange", "darkblue")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size = 12, color = colors),
        legend.position = "none") +
  ylab("Pearson's r") +
  xlab("") +
  scale_y_continuous(sec.axis = sec_axis(trans=~.*1, name="")) +
  geom_hline(yintercept = 0.184, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = -0.184, color = "grey", linetype = "dashed")

dir.create("./Figures/TCR/015_Linegraph_correlation_ES_TCR", showWarnings = FALSE)
pdf(paste0("./Figures/TCR/015_Linegraph_correlation_ES_TCR/Fig2a_", Source, "_linegraph_",
           "with_templates_and_TCR_clonality.pdf"), 
    width = 3.5, height = 4.5)
plot(plot)
dev.off()

