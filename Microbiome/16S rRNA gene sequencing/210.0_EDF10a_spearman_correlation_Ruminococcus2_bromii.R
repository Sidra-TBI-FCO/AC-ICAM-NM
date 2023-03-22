# common samples 

# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx"))

# set parameters

# load data
load("./Processed_Data/Survival Data/JSREP_NT_clinical_data.Rdata")
pcr = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/PCR_validation/Microbiome validation data all WGS cohort complete.xlsx")
samples.list = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/PCR_validation/list_of_samples_each_batch.xlsx")
load("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.0_Genus_full_abundance_filtered_10_percent_abundance_yes_0.01_T_246_samples.Rdata")

colSums(Genus_full_abundance)

Genus_full_abundance = t(Genus_full_abundance)
Genus_full_abundance = as.data.frame(Genus_full_abundance)

#Genus_full_abundance$`D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2`

Genus_full_abundance = Genus_full_abundance[which(rownames(Genus_full_abundance) %in% pcr$Sample),]

# subset PCR 
pcr = pcr[!is.na(pcr$`WGS.Species.=.Bromii`),]
pcr = pcr[!is.na(pcr$PCR),]

# 16s rum2
df.16s = pcr[,c(2,4,7)]
colnames(df.16s)[3] = "RA"
df.16s$platform = "16s"
df.16s$RA = NA
df.16s$RA = Genus_full_abundance$`D_0__Bacteria D_1__Firmicutes D_2__Clostridia D_3__Clostridiales D_4__Ruminococcaceae D_5__Ruminococcus 2`[match(pcr$Sample, rownames(Genus_full_abundance))]

# WGS 
df.wgs = pcr[,c(2,4,8)]
colnames(df.wgs)[3] = "RA"
df.wgs$platform = "WGS"

df.wgs$RA = df.wgs$RA/100

df = rbind(df.16s, df.wgs)

df$PCR = factor(df$PCR, levels = c("0" ,"1", "2", "3", "4"))

df$platform = factor(df$platform, levels = c("16s", "WGS"))


# correlation 

df.2 = cbind(df.16s, df.wgs)

colnames(df.2)[3] = "16s.Rum2"
colnames(df.2)[7] = "WGS.Bromii"

df.2 = df.2[,-c(4:6,8)]

# spearman correlation scatter blot 
svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/WGS_MB/spearman_correlation_RA_16s_rum2_wgs_bromii.svg"), width = 4, height = 4) 
plot = ggplot(df.2, aes(x= WGS.Bromii, y = `16s.Rum2`)) + 
  geom_point() +
  #scale_color_manual(values = c("ICR Low" = "blue", "ICR High" = "red")) +
  stat_cor(method = "spearman", size = 6) +
  #scale_y_log10() +
  #scale_x_log10() +
  theme_bw() +
  xlab("WGS bromii relative abundance") +
  ylab("16s Ruminococcus 2 relative abundance") +# TCR productive clonality
  theme(axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        aspect.ratio = 1/1) +
  geom_smooth(method="lm")

plot(plot)
dev.off()


dev.off()



