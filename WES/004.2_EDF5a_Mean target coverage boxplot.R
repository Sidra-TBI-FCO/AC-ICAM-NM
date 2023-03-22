
# Setup environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))

metrics = read.csv("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/QC Targeted Panel/table_csyl_29032020.csv", stringsAsFactors = FALSE)

#metrics = read.csv("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/QC Targeted Panel/table_slaj.csv", stringsAsFactors = FALSE)
#metrics = read.csv("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/OLD/QC Targeted Panel/table_fuqx.csv", stringsAsFactors = FALSE)

library(ggplot2)
library(ggrepel)
library(stringr)

metrics$tissue_type = str_sub(metrics$Sample, -1 , -1)
metrics$tissue_type[which(metrics$tissue_type %in% c("M", "1", "2"))] = "LM"
medians = aggregate(The.mean.coverage.of.targets. ~ tissue_type, metrics[,c("tissue_type", "The.mean.coverage.of.targets.")], FUN = median)
metrics$Sample_ID = gsub("WES_COAD_LUMC_SIDRA_", "", metrics$Sample)
#medians = medians$Mean.target.coverage
IQR = aggregate(The.mean.coverage.of.targets. ~ tissue_type, metrics[,c("tissue_type", "The.mean.coverage.of.targets.")], FUN = IQR)


pos = position_jitter(width = 0.5, seed = 1)
name = metrics$Sample_ID

metrics$tissue_type = factor(metrics$tissue_type, levels = c("N", "T", "LM"))
#metrics = metrics[-which(metrics$tissue_type == "LM"),]
levels(metrics$tissue_type) = c("Normal", "Tumor", "Tumor")

plot = ggplot(metrics, aes(x=tissue_type, y= The.mean.coverage.of.targets., fill=tissue_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.1) +
  scale_y_continuous(limits = c(0, 200)) +
  #geom_label_repel(label = name, position = pos) +
  scale_fill_manual(values = c("#E2F0D9", "#FBE5D6", "#FBE5D6")) +
  theme_bw() +
  ylab("Mean target coverage") +
  xlab("Tissue type") +
  theme(axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title = element_blank())
  #geom_text(data = medians, aes(x = tissue_type, y = The.mean.coverage.of.targets., label = The.mean.coverage.of.targets.), 
            #size = 3, vjust = -1.5)
  
dir.create("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/QC Targeted Panel/Boxplots", showWarnings = FALSE)
pdf("./WES_HPC/SDR_WH_WES_CGL_JSREP_DAVIDE/4_Post-Alignment_QC_WES/QC Targeted Panel/Boxplots/29032020_Mean_target_coverage_boxplot.pdf", width = 3.5, height = 3.5)
plot(plot)
dev.off()
