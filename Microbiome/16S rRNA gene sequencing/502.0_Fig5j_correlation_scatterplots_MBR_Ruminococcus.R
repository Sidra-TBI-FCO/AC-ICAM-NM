# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# load packages
required.packages = c("ggplot2", "reshape2", "ComplexHeatmap", "circlize", "dendsort", "ggpubr", "stringr")
ipak(required.packages)

# load data
#load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")
tumor.RA = Genus_full_abundance
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_N_246_samples_based_on_normal.Rdata")
normal.RA = Genus_full_abundance
rm(Genus_full_abundance)
load("./Processed_Data/Microbiome/External/9_August/ACICAM_Microbiome_Risk_scores.Rdata")
tumor.risk = Training
rm(Training, Validation)
normal.risk = read.csv("./Processed_Data/Microbiome/External/9_August/ACICAM_tumor_based_on_normal_without_5_genus_risk_Score.csv")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v3_clean_All_Immune_gene_signatures_table.Rdata")

# fix rownames
rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1,3)
immune_sig_df = immune_sig_df[which(rownames(immune_sig_df) %in% tumor.risk$Patient_ID),]

# risk tables
tumor.risk$samples = str_pad(tumor.risk$samples, 3, pad = 0)
normal.risk$samples = str_pad(normal.risk$samples, 3, pad = 0)

# add risk score
immune_sig_df$`tumor.risk` = tumor.risk$prediction[match(rownames(immune_sig_df), tumor.risk$samples)]
immune_sig_df$`normal.risk` = normal.risk$prediction[match(rownames(immune_sig_df), normal.risk$samples)]

# RA
colnames(tumor.RA) = substring(colnames(tumor.RA), 1, 3)
colnames(normal.RA) = substring(colnames(normal.RA), 1, 3)

tumor.RA = as.data.frame(t(tumor.RA))
normal.RA = as.data.frame(t(normal.RA))

# add ruminococcus2
immune_sig_df$`tumor.Rum2` = tumor.RA$`D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2`[match(rownames(immune_sig_df), rownames(tumor.RA))]
immune_sig_df$`normal.Rum2` = normal.RA$`D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2`[match(rownames(immune_sig_df), rownames(normal.RA))]


### scatter plot spearman  
sig1 = "tumor.Rum2" # x #tumor.Rum2
sig2 = "normal.Rum2" # y #tumor.risk,  normal.Rum2

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/scatter_plots/scatter_pearson_",sig1,"_",sig2,"_normal_risk_without_5_genera.pdf"), width = 3, height = 3) 
plot = ggplot(immune_sig_df, aes(x= immune_sig_df[,sig1], y = immune_sig_df[,sig2])) + 
  geom_point(size = 0.5) +
  stat_cor(method = "spearman", size = 5) +
  #scale_y_log10() +
  #scale_x_log10() +
  theme_bw() +
  xlab(paste0(sig1)) +
  ylab(paste0(sig2)) +# TCR productive clonality
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

plot(plot)
dev.off()

# exact value for immune_sig_df$tumor.Rum2, immune_sig_df$normal.Rum2
test = cor.test(immune_sig_df$tumor.Rum2, immune_sig_df$normal.Rum2, method = "spearman")
test$p.value

### scatter plot pearson 
sig1 = "normal.risk" # x 
sig2 = "tumor.risk" # y

test = cor.test(immune_sig_df$normal.risk, immune_sig_df$tumor.risk)
test$p.value

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/scatter_plots/scatter_pearson_",sig1,"_",sig2,"_normal_risk_without_5_genera.pdf"), width = 3, height = 3) 
plot = ggplot(immune_sig_df, aes(x= immune_sig_df[,sig1], y = immune_sig_df[,sig2])) + 
  geom_point(size = 0.5) +
  stat_cor(method = "pearson", size = 5) +
  #scale_y_log10() +
  #scale_x_log10() +
  theme_bw() +
  xlab("Normal risk score") +
  ylab("Tumor risk score") +# TCR productive clonality
  theme(axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

plot(plot)
dev.off()
