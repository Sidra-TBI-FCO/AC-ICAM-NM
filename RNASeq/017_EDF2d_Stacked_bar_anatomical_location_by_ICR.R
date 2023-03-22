
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
load("./Analysis/Trimmed_p/ICR Consensus Clustering/JSREP_ICR_cluster_assignment_k2-6.Rdata")
load("./Processed_Data/Survival Data/JSREP_clinical_data.Rdata")
excluded_df = read.csv("./Overview Sample Stocks/Meta_data/Excluded_patients.csv", stringsAsFactors = FALSE)
excluded_df$Patient_ID = str_pad(excluded_df$Patient_ID, 3, pad = "0")

clinical_data$tumour_anatomic_site = factor(clinical_data$tumour_anatomic_site, levels = c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                                           "colon transversum", "flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                                           "rectosigmoideum"))
clinical_data$primary_tumor_side = NA
clinical_data$primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("ceceum", "colon ascendens", "flexura hepatica", 
                                                                                 "colon transversum"))] = "Right sided"
clinical_data$primary_tumor_side[which(clinical_data$tumour_anatomic_site %in% c("flexura lienalis", "colon descendens", "colon sigmoideum",
                                                                                 "rectosigmoideum"))] = "Left sided"
clinical_data$primary_tumor_side = factor(clinical_data$primary_tumor_side, levels = c("Right sided", "Left sided"))

clinical_data$ICR = table_cluster_assignment$ICR_HML[match(clinical_data$Patient_ID, substring(rownames(table_cluster_assignment), 1, 3))]
clinical_data$ICR = factor(clinical_data$ICR, levels = c("ICR High", "ICR Medium", "ICR Low"))
clinical_data = clinical_data[-which(clinical_data$Patient_ID %in% excluded_df$Patient_ID),]
clinical_data = clinical_data[-which(is.na(clinical_data$ICR)),]

if(location == "All"){}else{
  clinical_data = clinical_data[which(clinical_data$primary_tumor_side == location),]
}

DF1 <- clinical_data %>%
  group_by(ICR) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF2 <- clinical_data %>%
  group_by(tumour_anatomic_site, ICR) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF3 <- clinical_data %>%
  group_by(tumour_anatomic_site, ICR) %>%
  summarise(count=n())

colors = c("ICR High" = "red",  "ICR Medium" = "green", 
           "ICR Low" = "blue")


plot2 = ggplot(DF2, aes(x = tumour_anatomic_site, y =perc*100, fill = ICR)) + geom_bar(stat="identity") +
  labs(x = "anatomic location", y = "Percentage", fill = "ICR", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")
        ) +
  scale_fill_manual(values= colors)


dir.create("./Figures/Trimmed_p/017.2_Stacked_barplot", showWarnings = FALSE)
pdf(paste0("./Figures/Trimmed_p/017.2_Stacked_barplot/017_Stacked_barchart_ICR_by_tumor_anatomic_location.pdf"),
    height = 6, width = 10)
plot(plot2)
dev.off()

plot = ggplot(DF1, aes(x = "", y =perc*100, fill = ICR)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "ICR", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black")) +
  scale_fill_manual(values= colors) + coord_polar("y", start=0)

#chisq_test(table(clinical_data$ICR, clinical_data$tumour_anatomic_site))

pdf(paste0("./Figures/Trimmed_p/017.2_Stacked_barplot/EDF2d_017.2_", location, "_348_ICR_Pie_Chart.pdf"),
    height = 7, width = 7)
plot(plot)
dev.off()

