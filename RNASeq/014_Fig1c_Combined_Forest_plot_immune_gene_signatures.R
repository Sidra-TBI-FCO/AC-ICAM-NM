
# Create forest plot of HR table generated in script 013

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("RColorBrewer", "forestplot", "stringr", "ggplot2")
ipak(required.packages)

# Set parameters
included_stages = "All"
Surv.cutoff.years = 20
exclude = c("Conpair_lower_90_percent", "non-epithelial")
Source = "ConsensusTME"

# Manual adjustment
order = c("ICR score", "Fibroblasts", "Endothelial", "Neutrophils", "Monocytes", "Macrophages", "Eosinophils", "Mast cells",
          "Cytotoxic cells", "NK cells", "T regulatory cells", "T cells gamma delta", "T cells CD8", "T cells CD4",
          "Immune Score", "Macrophages M1", "Dendritic cells", "Macrophages M2", "Plasma cells", "B cells")

# Load data
load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_", Source,"_OS_HR_table_exclude_",
            str_c(exclude, collapse = "_"), 
            "_stages_", str_c(included_stages, collapse = "_"), "_cutoff_", Surv.cutoff.years,".Rdata"))

HR_table$Signature = gsub("\\_", " ", HR_table$Signature)
HR_table$Signature[which(HR_table$Signature == "ICRscore")] = "ICR score"

HR_table = HR_table[order(HR_table$HR),]
rownames(HR_table) = HR_table$Signature
HR_table = HR_table[order,]
HR_table1 = HR_table

load(paste0("./Analysis/Trimmed_p/013_Immune_signatures_survival/013_", Source,"_DFS_Def1_HR_table_exclude_",
            str_c(exclude, collapse = "_"), 
            "_stages_", str_c(included_stages, collapse = "_"),
            "_cutoff_", Surv.cutoff.years,".Rdata"))

HR_table$Signature = gsub("\\_", " ", HR_table$Signature)
HR_table$Signature[which(HR_table$Signature == "ICRscore")] = "ICR score"

HR_table = HR_table[order(HR_table$HR),]
rownames(HR_table) = HR_table$Signature
HR_table = HR_table[order,]
HR_table2 = HR_table

HR_table = rbind(HR_table1, HR_table2)

## Forest plot seperate script
n_signatures = nrow(HR_table)
x = n_signatures + 2

HR.matrix = as.matrix(HR_table)
rownames(HR.matrix) = HR.matrix[,1]
HR.matrix = HR.matrix[,-c(1)]
mode(HR.matrix) = "numeric"

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR_table$HR[1:n_signatures]), NA),
    lower = c(NA,HR_table$CI_lower[c(1:n_signatures)], NA),
    upper = c(NA,HR_table$CI_upper[c(1:n_signatures)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")



HR_table$p_value = signif(HR_table$p_value, 3)
HR_table$HR = signif(HR_table$HR, 3)
tabletext<-cbind(
  c("Cell", as.character(HR_table$Signature)[c(1:n_signatures)]),
  c("p-value", HR_table$p_value[c(1:n_signatures)]),
  c("HR",      HR_table$HR[c(1:n_signatures)]))


dir.create("./Figures/Trimmed_p/014_Forest_plot_immune_signatures", showWarnings = FALSE)
pdf(file = paste0("./Figures/Trimmed_p/014_Forest_plot_immune_signatures/014_Forest_plot_", Source, "v4_immune_signatures_exclude_",
                  "Combined_OS_PFS_",
                  str_c(exclude, collapse = "_"),
                  "_cutoff_", Surv.cutoff.years, ".pdf"),
    height = 7, width = 4.5)
#height = 15, width = 7)
#height = 25, width = 15)
#height = 4, width = 3)

forestplot(mean = HR.matrix[,"HR"],
           lower = HR.matrix[,"CI_lower"],
           upper = HR.matrix[,"CI_upper"],
           labeltext = tabletext[-1,],
           new_page = FALSE,
           zero = 1,
           #is.summary=c(TRUE,rep(FALSE,n.cells),TRUE,rep(FALSE,n.cells),TRUE,FALSE),
           clip=c(0.001,55),
           xlog=TRUE,
           #xlim = c(0, 4),
           xlab = "HR",
           boxsize = .25,
           vertices = FALSE,
           col=fpColors(box="skyblue",line="skyblue"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 15),
                            ticks = gpar(fontsize = 15))
)
dev.off()


# Barplot 
HR_table$minus_log10_pval = -log10(HR_table$p_value)
HR_table$Survival_outcome = c(rep("OS", 20), rep("PFS", 20))
HR_table$Signature = factor(HR_table$Signature, levels = rev(order))
HR_table$Survival_outcome = factor(HR_table$Survival_outcome, c("PFS", "OS"))

plot = ggplot(HR_table, aes(x = Signature, y = minus_log10_pval, fill = Survival_outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values = c("OS" = "skyblue", "PFS" = "#FED8B1")) +
  coord_flip() + 
  scale_y_reverse() +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none") +
  xlab("") +
  ylab("-log10 p value") +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed")

pdf(paste0("./Figures/Trimmed_p/014_Forest_plot_immune_signatures/Combined_OS_PFS_Surv_cutoff_", Surv.cutoff.years,"years_minus_10_pval_barplot.pdf"),
    height = 4, width = 3)
plot(plot)
dev.off()
