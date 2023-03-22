
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "dplyr")
ipak(required.packages)

# Load data
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")

# Filter
clinical_data = clinical_data[which(clinical_data$Patient_ID %in% substring(colnames(RNASeq.QN.LOG2), 1, 3)),]

##### Sex

DF1 <- clinical_data %>%
  group_by(gender) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF1, aes(x = "", y =perc*100, fill = gender, label = count)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "gender", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, colour = "black")) +
  scale_fill_manual(values= c("FEMALE" = "#ff748c", "MALE" = "#74d2ff")) + coord_polar("y", start=0)

print(plot)

svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/AC_ICAM_clinical_gender.svg",  width = 3, height = 3)
print(plot)
dev.off()


################################################################################

# MSI

load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
MANTIS = MANTIS[which(MANTIS$Patient_ID %in% clinical_data$Patient_ID),]
MANTIS = MANTIS[-grep("LM", MANTIS$BAM.files),]

table(MANTIS$MSI)
MANTIS$MSI = factor(MANTIS$MSI, levels = c("MSS", "MSI-H"))


DF5 = MANTIS %>%
  group_by(MSI) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


plot5 = ggplot(DF5, aes(x = "", y =perc*100, fill = MSI, label = count)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "MSI", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.6)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("MSS" = "#FFFB01", "MSI-H" = "#A034F0")) + coord_polar("y", start=0)


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/AC_ICAM_clinical_MSI.svg", height = 3, width = 3)
plot(plot5)
dev.off()


###############################################################################

# stage

table(clinical_data$ajcc_pathologic_tumor_stage)

clinical = clinical_data

clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA"))] = "1"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"))] = "2"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = "3"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB"))] = "4"

table(clinical_data$ajcc_pathologic_tumor_stage)

DF3 = clinical_data %>%
  group_by(ajcc_pathologic_tumor_stage) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


plot3 = ggplot(DF3, aes(x = "", y =perc*100, fill = ajcc_pathologic_tumor_stage, label = count)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "ajcc_pathologic_tumor_stage", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.6)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("1" = "#55B78D", "2" = "#F3EB82",
                              "3" = "#F3C182", "4" = "#F66258")) + coord_polar("y", start=0)


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/AC_ICAM_clinical_stage.svg", height = 3, width = 3)
plot(plot3)
dev.off()

################################################################################

# age 

clinical = clinical_data

clinical$age_cat = NA
clinical$age_cat[which(clinical$age_at_initial_pathologic_diagnosis < 65 )] = "< 65"
clinical$age_cat[which(clinical$age_at_initial_pathologic_diagnosis >= 65 )] = ">= 65"

table(clinical$age_cat)

clinical$age_cat = factor(clinical$age_cat, levels = c("< 65", ">= 65"))

DF4 = clinical %>%
  group_by(age_cat) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


plot4 = ggplot(DF4, aes(x = "", y =perc*100, fill = age_cat, label = count)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "age_cat", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.6)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= c("< 65" = "#87CEEA", ">= 65" = "#90EE90")) + coord_polar("y", start=0)


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/AC_ICAM_clinical_age.svg", height = 3, width = 3)
plot(plot4)
dev.off()

################################################################################

# histology


clinical_data = clinical_data[which(clinical_data$Patient_ID %in% substring(colnames(RNASeq.QN.LOG2), 1, 3)),]

clinical_data$Tumor_morphology

clinical_data$histology = "Other"
clinical_data$histology[which(clinical_data$Tumor_morphology == "mucineus adenocarcinoom")] = "Mucineus"

table(clinical_data$histology)
clinical_data$histology = factor(clinical_data$histology, levels = c("Mucineus", "Other"))

DF5 <- clinical_data %>%
  group_by(histology) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF5, aes(x = "", y =perc*100, fill = histology, label = count)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "histology", face = "bold") +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, colour = "black")) +
  scale_fill_manual(values= c("Mucineus" = "#FF8C00", "Other" = "#AAA9AA")) + coord_polar("y", start=0)

print(plot)

svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/AC_ICAM_clinical_histology.svg",  width = 3, height = 3)
print(plot)
dev.off()











