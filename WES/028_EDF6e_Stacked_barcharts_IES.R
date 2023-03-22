
# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "dplyr"))

# Load data
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")
load("./Analysis/WES/002_Mutation_frequency/Nonsilent_mutation_frequency.Rdata")
load("./Processed_Data/WES/MANTIS/MANTIS.Rdata")
load("./Processed_Data/WES/MAF/finalMafFiltered_clean_primary_tumor_samples.Rdata")

POLE_mut_samples = finalMafFiltered$Patient_ID[which(finalMafFiltered$Hugo_Symbol == "POLE")]
  
# Prepare data

# Add mutational load
IES_df$Mut_load_cat = NA
IES_df$Mut_load = frequency_df$Nonsilent_mutational_burden_per_Mb[match(IES_df$Patient_ID,
                                                                   frequency_df$Patient_ID)]
IES_df$Mut_load_cat[which(IES_df$Mut_load >= 12)] = "Hypermutated"
IES_df$Mut_load_cat[which(IES_df$Mut_load < 12)] = "Non-hypermutated"

table(IES_df$Mut_load_cat)


#
IES_df$POLE = "WT"
IES_df$POLE[which(IES_df$Patient_ID %in% POLE_mut_samples)] = "MUT"

# Add MANTIS MSI
IES_df$MSI = MANTIS$MSI[match(IES_df$Patient_ID, MANTIS$Patient_ID)]

IES_df = IES_df[-which(is.na(IES_df$IES)),]


DF1 <- IES_df %>%
  group_by(IES, Mut_load_cat) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

DF2 <- IES_df %>%
  group_by(IES, MSI) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

Mut_cat_colors = c("Hypermutated" = "#EAAED0",
                   "Non-hypermutated" = "#A7EABD")

MSI_colors = c("MSI-H" = "purple",  "MSS" = "yellow")

plot = ggplot(DF1, aes(x = IES, y =perc*100, fill = Mut_load_cat)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values= Mut_cat_colors)

dir.create("./Figures/WES/028_Stacked_barchart_IES", showWarnings = FALSE)
svg("./Figures/WES/028_Stacked_barchart_IES/Mut_load_cat_IES_stacked_barchart.svg", 
    width = 5, height = 3, pointsize = 12)
plot(plot)
dev.off()

# MSI

plot = ggplot(DF2, aes(x = IES, y =perc*100, fill = MSI)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values= MSI_colors)

dir.create("./Figures/WES/028_Stacked_barchart_IES", showWarnings = FALSE)
svg("./Figures/WES/028_Stacked_barchart_IES/MSI_IES_stacked_barchart.svg", 
    width = 5, height = 3, pointsize = 12)
plot(plot)
dev.off()
