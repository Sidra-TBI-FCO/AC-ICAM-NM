
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("dplyr", "reshape2", "ggplot2", "survminer", "survival", "openxlsx", "stringr"))

#used = read.xlsx("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/PCR_validation/Microbiome validation data all WGS cohort complete.xlsx")
load("./Processed_Data/Survival Data/JSREP_NT_clinical_data.Rdata")
load("./Analysis/Trimmed_p/012_Immune_gene_signatures_table/012_v2_All_Immune_gene_signatures_table.Rdata")
load("./Processed_Data/Microbiome/From_WGS/100_WGS_167_matrices.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")

clinical_data = clinical_data[which(clinical_data$ID %in% colnames(Species_WGS)),]

Species_WGS = as.data.frame(t(Species_WGS))

clinical_data$Fusobacterium = Species_WGS$s__Fusobacterium_nucleatum[match(clinical_data$ID, rownames(Species_WGS))]
clinical_data$Fusobacterium = clinical_data$Fusobacterium/100

c1 = clinical_data[which(clinical_data$Fusobacterium == 0.0000000),]
c2 = clinical_data[-which(clinical_data$ID %in% c1$ID),]

c1$group = "Negative"

m.fuso = median(c2$Fusobacterium)

c2$group = NA

c2$group[which(c2$Fusobacterium < m.fuso)] = "Low"
c2$group[which(c2$Fusobacterium >= m.fuso)] = "High"

clinical_data = rbind(c1, c2)

table(clinical_data$group)

clinical_data$group = factor(clinical_data$group, levels = c("Negative", "Low", "High"))
table(clinical_data$group)

rownames(immune_sig_df) = substring(rownames(immune_sig_df), 1, 4)

clinical_data$Tcell = immune_sig_df$`Bindea |  T cells`[match(clinical_data$ID, rownames(immune_sig_df))]

clinical_data$MSI = MANTIS$MSI[match(clinical_data$ID, MANTIS$Sample_ID)]

table(clinical_data$MSI)
clinical_data$MSI = factor(clinical_data$MSI, levels = c("MSS", "MSI-H"))

clinical = clinical_data
clinical$Tcell.group = cut(clinical$Tcell, quantile(clinical$Tcell, prob = seq(0, 1, length = 4), type = 5))

clinical_data$tcell.group = NA

first = quantile(clinical_data$Tcell)[2]
second = quantile(clinical_data$Tcell)[3]
third = quantile(clinical_data$Tcell)[4]

clinical_data$tcell.group[which(clinical_data$Tcell < first)] = 1
clinical_data$tcell.group[which(clinical_data$Tcell >= first & clinical_data$Tcell < second)] = 2
clinical_data$tcell.group[which(clinical_data$Tcell >= second & clinical_data$Tcell < third)] = 3
clinical_data$tcell.group[which(clinical_data$Tcell >= third)] = 4

table(clinical_data$tcell.group)

# dir.create("./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/fusobacterium_analysis")
# save(clinical_data, file = "./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/fusobacterium_analysis/clinical_data_with_fusobacterium_grroups.Rdata")

# 
matrix = as.data.frame(t(RNASeq.QN.LOG2))
rownames(matrix) = substring(rownames(matrix), 1, 4)
clinical_data$CD163 = matrix$CD163[match(clinical_data$ID, rownames(matrix))]

clinical_data$group2 = NA
clinical_data$group2[which(clinical_data$group %in% c("Low", "High"))] = "Positive"
clinical_data$group2[which(clinical_data$group == "Negative")] = "Negative"

table(clinical_data$group2)

# quantile(clinical_data$Tcell)[2]
# max(clinical_data$Tcell)

my_comparison = list(c("Low", "High"), c("Negative", "Low"), c("Negative", "High"))

# boxplot 

pdf(file = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/boxplot/fusobacterium_nucleatum_tcell_two_groups_v2.pdf"),
    height = 3,width = 2)

ggplot(data = clinical_data, aes(x= group2, y= Tcell, fill = group2)) +
  geom_boxplot(outlier.colour = NA) +
  scale_fill_manual(values = c("#FFA07A", "#FF6347")) +
  scale_y_continuous("T cell") + # CD8 T cells
  #ggtitle(" ") +
  stat_compare_means(method = "t.test") +
  #ylim(2.5,20) +
  geom_point(position=position_jitterdodge(),alpha=1, size = 0.5) + 
  theme(plot.title = element_text(size=8, face = "bold")) + theme(axis.title = element_text(size = 12,  colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey", size = 0.2),
        panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
        axis.text.x = element_text(angle = 360, colour = "black", size = 10),
        axis.text.y = element_text(angle = 360, colour = "black", size = 10)) +
  theme(axis.line = element_line(color= "black", size = 0.4)) +
  guides(fill=FALSE)

dev.off()

###

clinical_data$tcell.group = factor(clinical_data$tcell.group, levels = c("1", "2", "3", "4"))

DF2 <- clinical_data %>%
  group_by(group, tcell.group) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

# DF2
pdf(file = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/boxplot/barplot_fuso_tcell_overall.pdf"),
    height = 3,width = 3)

ggplot(DF2, aes(x = group, y =perc*100, fill = tcell.group)) + geom_bar(stat="identity") +
  labs(x = "group", y = "Percentage", fill = "tcell.group", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 360,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = "black")
  ) +
  scale_fill_manual(values= c("#AFEEEE", "#48D1CC", "#5F9EA0", "#008B8B"))
dev.off()

class(clinical_data$group)
class(clinical_data$tcell.group)

table(clinical_data$group, clinical_data$tcell.group)
chisq.test(clinical_data$group, clinical_data$tcell.group)

# boxplot 

# svg(filename = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/boxplot/cd163_fusobacterium_groups_overall_v2.svg"),
#     height = 3,width = 2)
# 
# ggplot(data = clinical_data, aes(x= group, y= CD163, fill = group)) +
#   geom_boxplot(outlier.colour = NA) +
#   scale_fill_manual(values = c("#FFA07A", "#FF8C00", "#FF6347")) +
#   scale_y_continuous("CD163 gene") + 
#   #ggtitle(" ") +
#   stat_compare_means(comparisons = my_comparison, method = "t.test", size= 3) +
#   #ylim(2.5,20) +
#   geom_point(position=position_jitterdodge(),alpha=1, size = 0.5) + 
#   theme(plot.title = element_text(size=8, face = "bold")) + theme(axis.title = element_text(size = 12,  colour = "black")) + 
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(colour = "grey", size = 0.2),
#         panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
#         axis.text.x = element_text(angle = 360, colour = "black", size = 10),
#         axis.text.y = element_text(angle = 360, colour = "black", size = 10)) +
#   theme(axis.line = element_line(color= "black", size = 0.4)) +
#   guides(fill=FALSE)
# 
# dev.off()

# boxplot 

# svg(filename = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/boxplot/Tcells_fusobacterium_groups.svg"),
#     height = 3,width = 2)
# 
# ggplot(data = clinical_data, aes(x= group, y= Tcell, fill = group)) +
#   geom_boxplot(outlier.colour = NA) +
#   scale_fill_manual(values = c("#FFA07A", "#FF8C00", "#FF6347")) +
#   scale_y_continuous("T cells") + # CD8 T cells
#   #ggtitle(" ") +
#   #stat_compare_means(comparisons = my_comparison, method = "t.test",  label.x = 1.5, size = 3) +
#   #ylim(2.5,20) +
#   geom_point(position=position_jitterdodge(),alpha=1, size = 0.5) + 
#   theme(plot.title = element_text(size=8, face = "bold")) + theme(axis.title = element_text(size = 12,  colour = "black")) + 
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(colour = "grey", size = 0.2),
#         panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
#         axis.text.x = element_text(angle = 360, colour = "black", size = 10),
#         axis.text.y = element_text(angle = 360, colour = "black", size = 10)) +
#   theme(axis.line = element_line(color= "black", size = 0.4)) +
#   guides(fill=FALSE)
# 
# dev.off()


# boxplot 

# svg(filename = paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/boxplot/fusobacterium_nucleatum_cd163_two_groups.svg"),
#     height = 3,width = 2)
# 
# ggplot(data = clinical_data, aes(x= group2, y= CD163, fill = group2)) +
#   geom_boxplot(outlier.colour = NA) +
#   scale_fill_manual(values = c("#FFA07A", "#FF6347")) +
#   scale_y_continuous("CD163 genes") + # CD8 T cells
#   #ggtitle(" ") +
#   stat_compare_means(method = "t.test") +
#   #ylim(2.5,20) +
#   geom_point(position=position_jitterdodge(),alpha=1, size = 0.5) + 
#   theme(plot.title = element_text(size=8, face = "bold")) + theme(axis.title = element_text(size = 12,  colour = "black")) + 
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(colour = "grey", size = 0.2),
#         panel.grid.minor.y = element_line(colour = "grey", size = 0.2),
#         axis.text.x = element_text(angle = 360, colour = "black", size = 10),
#         axis.text.y = element_text(angle = 360, colour = "black", size = 10)) +
#   theme(axis.line = element_line(color= "black", size = 0.4)) +
#   guides(fill=FALSE)
# 
# dev.off()
