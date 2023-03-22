
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("ggplot2", "ggpubr", "dplyr")                                                                   
ibiopak(required.bioconductor.packages)

#load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")

colnames(Merged_dataset)

Merged_dataset$IES = IES_df$IES[match(Merged_dataset$Patient_ID, IES_df$Patient_ID)]

Merged_dataset$ajcc_pathologic_tumor_stage = factor(Merged_dataset$ajcc_pathologic_tumor_stage, levels = c("1", "2", "3", "4"))

DF1 <- Merged_dataset %>%
  group_by(IES, ajcc_pathologic_tumor_stage) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF1 = DF1[!is.na(DF1$IES),]

DF1$label = NA
DF1$label = DF1$perc*100
DF1$label = signif(DF1$perc*100, digits = 2)
DF1$label = paste0(DF1$label, "%")

Merged_dataset$tumor_stage = NA
Merged_dataset$tumor_stage[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% c(3, 4))] = "Stage III&IV"
Merged_dataset$tumor_stage[which(Merged_dataset$ajcc_pathologic_tumor_stage %in% c(1, 2))] = "Stage I&II"

Merged_dataset$IES_new = NA
Merged_dataset$IES_new[which(Merged_dataset$IES %in% c("IES1", "IES2"))] = "IES1+2"
Merged_dataset$IES_new[which(Merged_dataset$IES %in% c("IES3", "IES4"))] = "IES3+4"

table(Merged_dataset$IES_new, Merged_dataset$ajcc_pathologic_tumor_stage)
tbl = table(Merged_dataset$IES_new, Merged_dataset$ajcc_pathologic_tumor_stage)
tbl
chisq.test(tbl)
fisher.test(tbl)
chisq.test(tbl[c("IES1", "IES2"),])
fisher.test(tbl[c("IES1", "IES2"),])
chisq.test(tbl[c("IES3", "IES4"),])
fisher.test(tbl[c("IES3", "IES4"),])

tbl3 = table(Merged_dataset$IES_new, Merged_dataset$tumor_stage)
tbl3
chisq.test(tbl3)
fisher.test(tbl3)

Merged_dataset$tumor_stage_IV = as.character(Merged_dataset$ajcc_pathologic_tumor_stage)
Merged_dataset$tumor_stage_IV[which(Merged_dataset$tumor_stage_IV == "4")] = "Stage IV"
Merged_dataset$tumor_stage_IV[which(Merged_dataset$tumor_stage_IV %in% c("1", "2", "3"))] = "Stage I, II, III"

tbl4 = table(Merged_dataset$IES_new, Merged_dataset$tumor_stage_IV)
tbl4
chisq.test(tbl4)

Merged_dataset2 = Merged_dataset
tbl = table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage)
chisq.test(tbl)
fisher.test(tbl)
chisq.test(tbl[c("IES1", "IES2"),])
fisher.test(tbl[c("IES1", "IES2"),])
chisq.test(tbl[c("IES3", "IES4"),])
fisher.test(tbl[c("IES3", "IES4"),])

chisq.test(table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage))

Merged_dataset_new = Merged_dataset
Merged_dataset_new= Merged_dataset_new[-which(Merged_dataset_new$ICR_cluster == "ICR Medium"),]
Merged_dataset_new$Stage = as.character(Merged_dataset_new$ajcc_pathologic_tumor_stage)
Merged_dataset_new$Stage[which(Merged_dataset_new$Stage %in% c("1", "2", "3"))] = "Stage I, II, III"
Merged_dataset_new$Stage[which(Merged_dataset_new$Stage %in% c("4"))] = "Stage IV"

table(Merged_dataset_new$ICR_cluster, Merged_dataset_new$Stage)
fisher.test(table(Merged_dataset_new$ICR_cluster, Merged_dataset_new$Stage))
chisq.test(table(Merged_dataset_new$ICR_cluster, Merged_dataset_new$Stage))


table(Merged_dataset$IES, Merged_dataset$tumor_stage)
chisq.test(table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage))

#write.csv(DF1, file = "./Microbiome_EA/counts_stage_distribution_IES.csv")

#dir.create("./Analysis_EA/Figures/stacked_barplot")
png("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/stage_distribution_IES_barchart_count.png", res = 600, width = 4.2, height = 4, units = "in")

ggplot(DF1, aes(x = IES, y =perc*100, fill = ajcc_pathologic_tumor_stage, label = count)) + geom_bar(stat="identity") +
  labs(x = "IES", y = "Percentage", fill = "ajcc_pathologic_tumor_stage", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")
  ) +
  scale_fill_manual(values= c("#55B78E", "#F3EA81", "#F3C182", "#F66258"))
dev.off()

class(Merged_dataset$IES)

Merged_dataset$IES = factor(Merged_dataset$IES, levels = c("IES1", "IES2", "IES3", "IES4"))
Merged_dataset$ajcc_pathologic_tumor_stage = factor(Merged_dataset$ajcc_pathologic_tumor_stage, levels = c("1", "2", "3", "4"))

table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage)
chisq.test(table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage))

DF1$IES = factor(DF1$IES, levels = c("IES1", "IES2", "IES3", "IES4"))

table(DF1$IES, DF1$ajcc_pathologic_tumor_stage)
chisq.test(table(Merged_dataset$IES, Merged_dataset$ajcc_pathologic_tumor_stage))

table(Merged_dataset$ajcc_pathologic_tumor_stage[which(Merged_dataset$IES == "IES4")])
table(Merged_dataset$ICR_cluster[which(Merged_dataset$IES == "IES4")])

