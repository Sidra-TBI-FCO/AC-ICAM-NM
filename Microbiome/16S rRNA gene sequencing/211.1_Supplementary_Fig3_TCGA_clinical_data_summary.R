
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "dplyr")
ipak(required.packages)

# Load data
clinical_data = read.csv("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/External_Data/TCGA_CLINICAL_DATA_CELL_2018_S1.csv", stringsAsFactors = FALSE)
#load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
load("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/External_Data/From_TCGA_Masterfile/Biolinks_MANTIS.Rdata")


colnames(clinical_data)


# Filter
clinical_data = clinical_data[which(clinical_data$bcr_patient_barcode %in% MANTIS$Patient_ID),]
MANTIS = MANTIS[which(MANTIS$Patient_ID %in% clinical_data$bcr_patient_barcode),]

##############
####### gender 

DF1 <- clinical_data %>%
  group_by(gender) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

#dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA")

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

svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA/clinical_gender.svg",  width = 3, height = 3)
print(plot)
dev.off()

################################################################################
################################################################################
# tumor stage

table(clinical_data$ajcc_pathologic_tumor_stage)

clinical = clinical_data

clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA"))] = "1"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"))] = "2"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = "3"
clinical$ajcc_pathologic_tumor_stage[which(clinical$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB"))] = "4"

table(clinical$ajcc_pathologic_tumor_stage)

DF3 = clinical %>%
  group_by(ajcc_pathologic_tumor_stage) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF3= DF3[-c(1,2),]

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


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA/clinical_stage.svg", height = 3, width = 3)
plot(plot3)
dev.off()

################################################################################
# age

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


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA/clinical_age.svg", height = 3, width = 3)
plot(plot4)
dev.off()

################################################################################

# MSI

DF5 = MANTIS %>%
  group_by(MSI) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))


DF5 = DF5[-3,]


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


svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA/clinical_MSI.svg", height = 3, width = 3)
plot(plot5)
dev.off()

################################################################################

# histology

clinical_tcga = read.csv("../NGS_Data_TCGA_COAD_Jessica/Processed_Data/External_Data/clinical_data_from_PancancerAtlas.csv", stringsAsFactors = FALSE)

clinical_tcga = clinical_tcga[which(clinical_tcga$submitter_id %in% clinical_data$bcr_patient_barcode),]
clinical_tcga = clinical_tcga[!duplicated(clinical_tcga$submitter_id),]

clinical_tcga$primary_diagnosis

clinical_tcga$histology = "Other"
clinical_tcga$histology[which(clinical_tcga$primary_diagnosis == "Mucinous adenocarcinoma")] = "Mucineus"

table(clinical_tcga$histology)
clinical_tcga$histology = factor(clinical_tcga$histology, levels = c("Mucineus", "Other"))

DF7 <- clinical_tcga %>%
  group_by(histology) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF7, aes(x = "", y =perc*100, fill = histology, label = count)) + geom_bar(stat="identity") +
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

svg("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/TCGA/clinical_histology.svg",  width = 3, height = 3)
print(plot)
dev.off()



