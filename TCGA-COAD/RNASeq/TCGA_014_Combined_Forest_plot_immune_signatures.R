

# Create forest plot of HR table generated in script 013

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages <- c("RColorBrewer", "forestplot", "stringr", "ggplot2")
ipak(required.packages)

# Set parameters
Surv.cutoff.years = 20
geneset = "ConsensusTME_COAD"

# Manual adjustment
order = c("ICR score", "Fibroblasts", "Endothelial", "Neutrophils", "Monocytes", "Macrophages", "Eosinophils", "Mast cells",
          "Cytotoxic cells", "NK cells", "T regulatory cells", "T cells gamma delta", "T cells CD8", "T cells CD4",
          "Immune Score", "Macrophages M1", "Dendritic cells", "Macrophages M2", "Plasma cells", "B cells")

# Load data
load(paste0("./Analysis/013_Immune_signatures_survival/013_", geneset, "_OS_HR_table_",
     "cutoff_", Surv.cutoff.years,".Rdata"))

HR_table$Signature = gsub("_", " ", HR_table$Signature)

HR_table = HR_table[order(HR_table$HR),]
rownames(HR_table) = HR_table$Signature
HR_table = HR_table[order,]
HR_table1 = HR_table

load(paste0("./Analysis/013_Immune_signatures_survival/013_", geneset, "_PFI_HR_table_",
            "cutoff_", Surv.cutoff.years,".Rdata"))

HR_table$Signature = gsub("_", " ", HR_table$Signature)

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


dir.create("./Figures/014_Forest_plot_immune_signatures", showWarnings = FALSE)
pdf(file = paste0("./Figures/014_Forest_plot_immune_signatures/014_Forest_plot_", geneset, "_immune_signatures_",
                  "Combined_OS_PFI",
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
           col=fpColors(box="#FED8B1",line="#FED8B1"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 15),
                            ticks = gpar(fontsize = 15))
)
dev.off()

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

png(paste0("./Figures/014_Forest_plot_immune_signatures/Combined_OS_PFI_Surv_cutoff_", Surv.cutoff.years,"years_minus_10_pval_barplot.png"), res = 600,
    units = "in", height = 4, width = 3)
plot(plot)
dev.off()
