
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data_TCGA_COAD_Jessica"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "dplyr")
ipak(required.packages)

# Set parameters
location = "All" # "All" or "Right sided" or "Left sided"

# Load data
load("./Analysis/016_CMS_Classification/Biolinks_Rfcms.Rdata")
load("./Processed_Data/RNASeq/004_TCGA_COAD_RNASeq_GDC_Biolinks_Normalized_TP_Extracted_Whitelist_Filtered.Rdata")
clinical_data = read.csv("./Processed_Data/External_Data/clinical_data_from_PancancerAtlas.csv", stringsAsFactors = FALSE)

clinical_data$site_of_resection_or_biopsy = factor(clinical_data$site_of_resection_or_biopsy, levels = c("Cecum", "Ascending colon", "Hepatic flexure of colon", 
                                                                                           "Transverse colon", "Splenic flexure of colon", "Descending colon", 
                                                                                           "Sigmoid colon", "Rectosigmoid junction", "Colon, NOS"))
clinical_data$primary_tumor_side = NA
clinical_data$primary_tumor_side[which(clinical_data$site_of_resection_or_biopsy %in% c("Cecum", "Ascending colon", "Hepatic flexure of colon", 
                                                                                 "Transverse colon"))] = "Right sided"
clinical_data$primary_tumor_side[which(clinical_data$site_of_resection_or_biopsy %in% c("Splenic flexure of colon", "Descending colon", 
                                                                                 "Sigmoid colon", "Rectosigmoid junction"))] = "Left sided"
clinical_data$primary_tumor_side = factor(clinical_data$primary_tumor_side, levels = c("Right sided", "Left sided"))

Rfcms$Patient_ID = substring(rownames(Rfcms), 1, 12)
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed or indeterminate"
Rfcms$anatomic_location = NA

Rfcms$anatomic_location = clinical_data$site_of_resection_or_biopsy[match(Rfcms$Patient_ID,
                                                                   clinical_data$submitter_id)]
Rfcms$primary_tumor_side = clinical_data$primary_tumor_side[match(Rfcms$Patient_ID,
                                                                  clinical_data$submitter_id)]

Rfcms = Rfcms[-which(is.na(Rfcms$anatomic_location)),]

if(location == "All"){}else{
  Rfcms = Rfcms[which(Rfcms$primary_tumor_side == location),]
}

Rfcms = Rfcms[which(rownames(Rfcms) %in% colnames(filtered.norm.RNAseqData)),]

DF1 <- Rfcms %>%
  group_by(RF.predictedCMS) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF2 <- Rfcms %>%
  group_by(anatomic_location, RF.predictedCMS) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF3 <- Rfcms %>%
  group_by(anatomic_location, RF.predictedCMS) %>%
  summarise(count=n())

colors = c("CMS1" = "#FF9F21",  "CMS2" = "#0074AF", 
           "CMS3" = "#E97AA8", "CMS4" = "#009E74",
           "mixed or indeterminate" = "lightgrey")

plot = ggplot(DF1, aes(x = "", y =perc*100, fill = RF.predictedCMS)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "CMS", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors) + coord_polar("y", start=0)

plot2 = ggplot(DF2, aes(x = anatomic_location, y =perc*100, fill = RF.predictedCMS)) + geom_bar(stat="identity") +
  labs(x = "anatomic location", y = "Percentage", fill = "CMS", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors)

plot3 = ggplot(DF3, aes(x = anatomic_location, y =count, fill = RF.predictedCMS)) + geom_bar(stat="identity") +
  labs(x = "anatomic location", y = "Count", fill = "CMS", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors)

dir.create("./Figures/017_CMS_Circle_plot", showWarnings = FALSE)
png(paste0("./Figures/017_CMS_Circle_plot/017_", location, "_Biolinks_CMS_Pie_Chart.png"),
    height = 7, width = 7, units = "in", res = 600)
plot(plot)
dev.off()

png(paste0("./Figures/017_CMS_Circle_plot/017_Stacked_barchart_by_tumor_anatomic_location.png"),
    height = 6, width = 10, units = "in", res = 600)
plot(plot2)
dev.off()

png(paste0("./Figures/017_CMS_Circle_plot/017_Stacked_barchart_count_by_tumor_anatomic_location.png"),
    height = 6, width = 10, units = "in", res = 600)
plot(plot3)
dev.off()
