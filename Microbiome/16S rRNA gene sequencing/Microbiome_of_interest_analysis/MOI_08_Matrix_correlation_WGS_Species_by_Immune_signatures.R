
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))


# Set parameters
Species = "s__Fusobacterium_nucleatum"

# Load data
load(paste0("./Analysis/Microbiome/Microbiome_Of_Interest_Analysis/004_correlation_with_gene_signatures/immune_sig_df_",
            Species, ".Rdata"))

colnames(immune_sig_df)[104] = Species

# Plot

# Correlation panel
reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x), col = "blue")
  r <- round(cor(x, y, method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.2, 0.9, txt)
}



dir.create("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_08_Matrix_Correlation_Immune_Trait_and_Microbiome", showWarnings = FALSE)
pdf(file = paste0("./Figures/Microbiome/Microbiome_Of_Interest_Plots/MOI_08_Matrix_Correlation_Immune_Trait_and_Microbiome/",
                  Species, "_by_2_signatures.pdf"), width = 8, height = 8)
# Create the plots
pairs(immune_sig_df[,c(Species, "Expression Signature - TREM1 data", "Expression Signature - CD68", "Expression Signature - IL8 21978456")], 
      lower.panel = NULL,
      upper.panel = reg)

dev.off()
