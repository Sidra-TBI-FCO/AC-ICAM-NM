
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr", "ggplot2", "ggpubr", "dplyr")
ipak(required.packages)

# Set parameters
location = "Left sided" # "All" or "Right sided" or "Left sided"
exclude = c("Conpair_lower_90_percent", "non-epithelial")

# Load data
load("./Analysis/Trimmed_p/016_CMS_Classification/Rfcms.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

clinical_data$tumour_anatomic_site = factor(clinical_data$tumour_anatomic_site, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                                                         "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",                                                                                                         "rectosigmoideum"))
clinical_data$primary_tumor_side = NA
clinical_data$primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                                       "colon transversum"))] = "Right sided"
clinical_data$primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                                       "rectosigmoideum"))] = "Left sided"
clinical_data$primary_tumor_side = factor(clinical_data$primary_tumor_side, levels = c("Right sided", "Left sided"))

Rfcms$Patient_ID = substring(rownames(Rfcms), 1, 3)
Rfcms$RF.predictedCMS[which(is.na(Rfcms$RF.predictedCMS))] = "mixed or indeterminate"
Rfcms$anatomic_location = NA

Rfcms$anatomic_location = clinical_data$tumour_anatomic_site[match(Rfcms$Patient_ID,
                                                                   clinical_data$Patient_ID)]
Rfcms$primary_tumor_side = clinical_data$primary_tumor_side[match(Rfcms$Patient_ID,
                                                                  clinical_data$Patient_ID)]

if(location == "All"){}else{
  Rfcms = Rfcms[which(Rfcms$primary_tumor_side == location),]
}
if(exclude == ""){
}else{
  Rfcms = Rfcms[-which(Rfcms$Patient_ID %in% excluded_df$Patient_ID[which(excluded_df$Reason.excluded %in% exclude)]), ]
}

DF1 <- Rfcms %>%
  group_by(RF.predictedCMS) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF2 <- Rfcms %>%
  group_by(anatomic_location, RF.predictedCMS) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

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

dir.create("./Figures/Trimmed_p/017_CMS_Circle_plot", showWarnings = FALSE)
pdf(paste0("./Figures/Trimmed_p/017_CMS_Circle_plot/017_EDF2d_", location, "_348_CMS_Pie_Chart.pdf"),
    height = 7, width = 7)
plot(plot)
dev.off()

pdf(paste0("./Figures/Trimmed_p/017_CMS_Circle_plot/017EDF2d_CMS_Stacked_barchart_by_tumor_anatomic_location.pdf"),
    height = 6, width = 10)
plot(plot2)
dev.off()

