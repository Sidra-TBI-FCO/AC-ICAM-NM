

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
source(paste0(toolbox.path, "/R scripts/ggkm_Jessica_Pancancer.R"))

ipak(c("stringr", "ggplot2", "ggpubr", "dplyr", "scales"))

# Set parameters
exclude_medium = "include_medium_cat"  # "include_medium_cat" "include_medium" or "exclude_medium"
Group.of.interest = "Immunoedited_ICR_cluster" # "Mutload_ICR_cluster" or "Immunoedited_ICR_cluster"
Surv.cutoff.years = 20
x = 2
Stages = "all stages" #"all stages" # c(3, 4) # "all stages"  or c(1, 2) or c(3, 4)
Stage_name = "All" #"Stage III&IV" #"All"

# Load data
load(paste0("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset_",
            exclude_medium, ".Rdata"))
MiXCR_orig = read.csv("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Diversity.csv", 
                      stringsAsFactors = FALSE)

# Prepare MiXCR data
MiXCR = MiXCR_orig
colnames(MiXCR)[1] = "Sample_ID"
MiXCR$Patient_ID = substring(MiXCR$Sample_ID, 1, 3)
MiXCR$Tissue = substring(MiXCR$Sample_ID, 4, 6)

MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),] # Filter only tumor samples
MiXCR = MiXCR[which(substring(MiXCR$Sample_ID, 7, 10) == "_TRB"),]
MiXCR = MiXCR[-which(is.na(MiXCR$Pielou)),]
MiXCR$TCR_Clonality = 1 - MiXCR$Pielou

MiXCR = MiXCR[-which(MiXCR$TCR_Clonality == -Inf),]

# Overview subgroups
Merged_dataset$TCR_MiXCR_clonality = MiXCR$TCR_Clonality[match(Merged_dataset$Patient_ID,
                                                               MiXCR$Patient_ID)]
table(Merged_dataset$Immunoedited)
table(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)
Merged_dataset$Immunoedited_ICR_cluster = paste(Merged_dataset$Immunoedited, Merged_dataset$ICR_cluster)
Merged_dataset$Mutload_ICR_cluster = paste(Merged_dataset$Mutation_cat, Merged_dataset$ICR_cluster)

table(Merged_dataset$ajcc_pathologic_tumor_stage)
if(Stages == "all stages"){
}else{
  Merged_dataset = Merged_dataset[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% Stages),]
}
table(Merged_dataset$ajcc_pathologic_tumor_stage, Merged_dataset$Immunoedited_ICR_cluster)

table(Merged_dataset$Immunoedited_ICR_cluster)

if(exclude_medium == "include_medium_cat"){
  levels = rev(c("immunoedited ICR High", "less immunoedited ICR High",
              "immunoedited ICR Medium", "less immunoedited ICR Medium",
              "immunoedited ICR Low", "less immunoedited ICR Low"))
}else{
  levels = rev(c("immunoedited ICR High", "less immunoedited ICR High",
             "immunoedited ICR Low", "less immunoedited ICR Low"))
}

Merged_dataset$Immunoedited_ICR_cluster = factor(Merged_dataset$Immunoedited_ICR_cluster,
                                                 levels = levels)

Merged_dataset$Mutload_ICR_cluster = factor(Merged_dataset$Mutload_ICR_cluster,
                                            levels = c("nonhypermutated ICR Low", "hypermutated ICR Low",
                                                       "nonhypermutated ICR High", "hypermutated ICR High"))

if(Group.of.interest == "Immunoedited_ICR_cluster" & exclude_medium %in% c("include_medium", "exclude_medium")){
  col_values = c("immunoedited ICR High" = "#FF3806",
                 "less immunoedited ICR High" = "#FF9F00",
                 "immunoedited ICR Low" = "#009CFC",
                 "less immunoedited ICR Low" = "#5233FC")
}
if(Group.of.interest == "Immunoedited_ICR_cluster" & exclude_medium %in% c("include_medium_cat")){
  col_values = c("immunoedited ICR High" = "#FF3806",
                 "less immunoedited ICR High" = "#FF9F00",
                 "immunoedited ICR Medium" = "darkgreen",
                 "less immunoedited ICR Medium" = "lightgreen",
                 "immunoedited ICR Low" = "#009CFC",
                 "less immunoedited ICR Low" = "#5233FC")
}

if(Group.of.interest == "Mutload_ICR_cluster"){
  col_values = c("nonhypermutated ICR Low" = "#84D976",
                 "hypermutated ICR Low" = "#7C7AED",
                 "nonhypermutated ICR High" = "#00871F",
                 "hypermutated ICR High" = "#8B16DD")
}


table(Merged_dataset$ICR_cluster)
nrow(Merged_dataset) - sum(is.na(Merged_dataset$TCR_productive_clonality))

nrow(Merged_dataset) - sum(is.na(Merged_dataset$TCR_MiXCR_clonality))
MiXCR_subset = Merged_dataset[-which(is.na(Merged_dataset$TCR_MiXCR_clonality)),]
table(MiXCR_subset$Immunoedited_ICR_cluster)

violin_plot = ggplot(Merged_dataset, aes(x = get(Group.of.interest), y = TCR_productive_clonality)) +
                                         #label = round(ICRscore, 1))) +
  geom_violin(aes(fill = get(Group.of.interest))) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = get(Group.of.interest))) +
  #geom_jitter(width = 0.1) + geom_text_repel() + stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
   #                                                                                    c("immunoedited ICR Low", "less immunoedited ICR Low")), method = "t.test") + 
  #ylim(0, 0.4) + # extra edit, not in main figure
  stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
                                        c("immunoedited ICR Medium", "less immunoedited ICR Medium"),
                                        c("immunoedited ICR Low", "less immunoedited ICR Low")),
                     method = "t.test") +
  theme_bw() +
  xlab("") +
  ylab("TCR productive clonality \n (ImmunoSeq)") +
  scale_fill_manual(values = col_values) +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 15, angle = 90,
         #                         vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.position = "none")

pdf(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/Fig4e_Violinplot_", Group.of.interest, "_", exclude_medium,"TCR_ImmunoSeq_ICR_by_category.pdf"),
    width = 5, height = 4)
plot(violin_plot)
dev.off()

test = t.test(Merged_dataset$TCR_productive_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "less immunoedited ICR Low")],
              Merged_dataset$TCR_productive_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "immunoedited ICR Low")])
test$p.value

test = t.test(Merged_dataset$TCR_productive_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "less immunoedited ICR High")],
              Merged_dataset$TCR_productive_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "immunoedited ICR High")])
test$p.value

violin_plot2 = ggplot(Merged_dataset, aes(x = get(Group.of.interest), y = TCR_MiXCR_clonality)) +
                                          #label = round(ICRscore, 1))) +
  geom_violin(aes(fill = get(Group.of.interest))) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = get(Group.of.interest))) +
  theme_bw() +
  stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
                                        c("immunoedited ICR Medium", "less immunoedited ICR Medium"),
                                        c("immunoedited ICR Low", "less immunoedited ICR Low")),
                     method = "t.test") +
  #geom_jitter(width = 0.1) + geom_text_repel() + stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
   #                                                                                    c("immunoedited ICR Low", "less immunoedited ICR Low")), method = "t.test") + 
  #ylim(0, 0.35) + # extra edit, not in main figure
  xlab("") +
  ylab("TCRB clonality \n (MiXCR)") +
  scale_fill_manual(values = col_values) +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 15, angle = 90,
        #                          vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.position = "none")

pdf(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/Fig4e_Violinplot_",Group.of.interest, "_", exclude_medium,"TCR_MiXCR_ICR_by_category.pdf"),
    width = 5, height = 4)
plot(violin_plot2)
dev.off()

test = t.test(Merged_dataset$TCR_MiXCR_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "less immunoedited ICR Low")],
              Merged_dataset$TCR_MiXCR_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "immunoedited ICR Low")])
test$p.value

test = t.test(Merged_dataset$TCR_MiXCR_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "less immunoedited ICR High")],
              Merged_dataset$TCR_MiXCR_clonality[which(Merged_dataset$Immunoedited_ICR_cluster == "immunoedited ICR High")])
test$p.value

medium_only = Merged_dataset[which(Merged_dataset$ICR_cluster == "ICR Medium"),]

violin_plot3 = ggplot(medium_only, aes(x = get(Group.of.interest), y = TCR_MiXCR_clonality)) +
  geom_violin(aes(fill = get(Group.of.interest))) +
  geom_boxplot(width=.1, outlier.shape = NA, aes(fill = get(Group.of.interest))) +
  theme_bw() +
  xlab("") +
  ylab("TCRB clonality \n (MiXCR)") +
  scale_fill_manual(values = col_values) +
  #stat_compare_means(method = "t.test") +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 15, angle = 90,
        #                          vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.position = "none")

png(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/EDF6i_Violinplot_ICR_Medium_only_",Group.of.interest, "_", exclude_medium,"TCR_MiXCR_ICR_by_category.png"),
    width = 2.3, height = 3, units = "in", res = 600)
plot(violin_plot3)
dev.off()

pdf(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/EDF6i_Violinplot_ICR_Medium_only_",Group.of.interest, "_", exclude_medium,"TCR_MiXCR_ICR_by_category.pdf"),
    width = 2.3, height = 3)
plot(violin_plot3)
dev.off()

plot4 = ggplot(Merged_dataset, aes(x = get(Group.of.interest), y = TCR_MiXCR_clonality)) +
  geom_boxplot(outlier.shape = NA, aes(fill = get(Group.of.interest))) +
  geom_jitter(size = 0.2, width = 0.2) +
  theme_bw() +
  xlab("") +
  ylab("TCRB clonality \n (MiXCR)") +
  #geom_smooth(method = "lm", se=TRUE, aes(group=1)) +
  scale_fill_manual(values = col_values) + 
  #stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
   #                                     c("immunoedited ICR Low", "less immunoedited ICR Low")),
    #                 method = "t.test") +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 15, angle = 90,
         #                          vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.position = "none")

png(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/Boxplot_", exclude_medium,"TCR_MiXCR_ICR_by_category.png"),
    width = 3, height = 3, units = "in", res = 600)
plot(plot4)
dev.off()


Merged_dataset$Immunoedited_ICR_cluster_num = as.numeric(Merged_dataset$Immunoedited_ICR_cluster)

cor.test(x = Merged_dataset$Immunoedited_ICR_cluster_num, y= Merged_dataset$TCR_MiXCR_clonality, method = "spearman")
cor.test(x = Merged_dataset$Immunoedited_ICR_cluster_num, y= Merged_dataset$TCR_productive_clonality, method = "spearman")

plot5 = ggplot(Merged_dataset, aes(x = get(Group.of.interest), y = Nonsilent_mutational_burden_per_Mb)) +
  geom_boxplot(outlier.shape = NA, aes(fill = get(Group.of.interest))) +
  geom_jitter(size = 0.2, width = 0.2) +
  theme_bw() +
  xlab("") +
  ylab("Nonsynonymous Mutational load \n (per Mb)") +
  scale_fill_manual(values = col_values) + 
  stat_compare_means(comparisons = list(c("immunoedited ICR High", "less immunoedited ICR High"),
                                        c("immunoedited ICR Low", "less immunoedited ICR Low")),
                     method = "t.test") +
  scale_y_log10(labels = format_format(scientific = FALSE)) +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black", size = 15, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15))

png(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/Boxplot_", exclude_medium,"Mutational_load_by_category.png"),
    width = 4.5, height = 6, units = "in", res = 600)
plot(plot5)
dev.off()


DF <- Merged_dataset %>%
  group_by(ajcc_pathologic_tumor_stage, get(Group.of.interest)) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot6 = ggplot(DF, aes(x = ajcc_pathologic_tumor_stage, y =perc*100, fill = get(Group.of.interest))) + geom_bar(stat="identity") +
  labs(x = "AJCC pathological stage", y = "Percentage", fill = "Immunoedited_ICR_cluster", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"), 
        legend.position = "none") +
  scale_fill_manual(values= col_values)

png(paste0("./Figures/WES/022f_Boxplots_ICRscore_IE_score/Stacked_barchart_", exclude_medium,"_Pathological_stage_by_category.png"),
    width = 4, height = 4, units = "in", res = 600)
plot(plot6)
dev.off()


## Save IES categories
IES_df = Merged_dataset
IES_df = IES_df[, c("Patient_ID", "Immunoedited_ICR_cluster", "MSI_status")]
IES_df$IES = NA
IES_df$IES[which(IES_df$Immunoedited_ICR_cluster == "less immunoedited ICR Low")] = "IES1"
IES_df$IES[which(IES_df$Immunoedited_ICR_cluster == "immunoedited ICR Low")] = "IES2"
IES_df$IES[which(IES_df$Immunoedited_ICR_cluster == "less immunoedited ICR High")] = "IES3"
IES_df$IES[which(IES_df$Immunoedited_ICR_cluster == "immunoedited ICR High")] = "IES4"

table(IES_df$MSI_status, IES_df$IES)
IES_df$MSI_status = NULL

save(IES_df, file = "./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")

