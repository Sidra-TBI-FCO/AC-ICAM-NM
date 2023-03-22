

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2"))

# Set parameters
per = "no" # 10 # 5 no
abundance = "no" # no  zyes_0.01 
Type = "Relative"
Mode = "all_samples" # "all_tumor" or "paired"
Tissue = "T"

# Load data
load(paste0("./Processed_Data/Microbiome/001_data_preparation/OTU_tables/", Type,"_abundancies/clean_", Mode, "_rank_level_abundancies.Rdata"))
poore_contaminants_list = read.csv('./Processed_Data/External/Microbiome/Combined_TableS6+S7_Poore.csv',
                                   stringsAsFactors = FALSE)

load("./Overview Sample Stocks/Overview_available_data_per_patient/Overview_26_April_2022.Rdata")
Cohort_246 = paste(Overview$Patient[which(Overview$Microbiome == "yes")], "", sep = "")

# Annotate
dim(Genus_full_abundance)

full = rownames(Genus_full_abundance)

df = data.frame(Full = full,
                Kingdom = gsub(".*\\D_0__", "", gsub("\\D_1__.*", "", full)),
                Phylum = gsub(".*\\D_1__", "", gsub("\\D_2__.*", "", full)),
                Class = gsub(".*\\D_2__", "", gsub("\\D_3__.*", "", full)),
                Order = gsub(".*\\D_3__", "", gsub("\\D_4__.*", "", full)),
                Genus = gsub(".*\\D_4__", "", gsub("\\D_5__.*", "", full)),
                Species = gsub(".*\\D_5__", "", gsub("\\D_6__.*", "", full)))

S6_contamin = poore_contaminants_list$Genera[which(poore_contaminants_list$Source == "S6")]
S7_contamin = poore_contaminants_list$Genera[which(poore_contaminants_list$Source == "S7")]

df$Contaminant = NA
df$Contaminant[which(df$Species %in% S6_contamin)] = "S6 Contaminant"
df$Contaminant[which(df$Species %in% S7_contamin)] = "S7 Contaminant"
  
table(df$Contaminant, exclude = NULL)

df_contaminants = df[-which(is.na(df$Contaminant)),]

df_contaminants$Category = poore_contaminants_list$Category[match(df_contaminants$Species,
                                                                  poore_contaminants_list$Genera)]

table(df_contaminants$Contaminant, df_contaminants$Category)

dir.create("./Analysis/Exploration_reviewers/Microbiome/026_Decontamination_prepare_overview", showWarnings = FALSE)
write.csv(df_contaminants, file = paste0("./Analysis/Exploration_reviewers/Microbiome/026_Decontamination_prepare_overview/", 
                                         "Potential_contamination_overview_AC_ICAM_", Type, "_abundancies_", Mode, "_percentage_", 
                                         per, "_", abundance ,".csv"),
          row.names = FALSE)


# Plotting the putative contaminants
contaminants = df_contaminants[which(df_contaminants$Category %in% c("MIXED EVIDENCE", "LIKELY CONTAMINANT") &
                                          df_contaminants$Contaminant == "S7 Contaminant"),]

contam = contaminants$Full

save(contam, file = "./Analysis/Exploration_reviewers/Microbiome/tumor_only_246/01.26_Decontamination_prepare_overview/Contam_35_to_exclude.Rdata")
## Plot relative abundance per sample for all contaminants

# Stacked barchart
Genus_full_abundance = Genus_full_abundance[c(contam, 
                                              rownames(Genus_full_abundance)[-which(rownames(Genus_full_abundance) %in% contam)]),]
Genus_full_abundance$Genus = rownames(Genus_full_abundance)
DF = melt(Genus_full_abundance, id = "Genus")
DF$Tissue = substring(DF$variable, 4,4)

DF = DF[which(DF$Tissue == Tissue),]

colors_df = data.frame(Genus = rownames(Genus_full_abundance), Color = "grey")

colors_df$Color[which(colors_df$Genus %in% contam)] = "red"

Contam_colors = colors_df$Color
names(Contam_colors) = colors_df$Genus

plot = ggplot(DF, aes(x = variable, y = value*100, fill = Genus)) + geom_bar(stat="identity") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "darkgrey"),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, colour = "black", size = 7),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none",
        axis.ticks = element_blank()) +
  scale_fill_manual(values = Contam_colors) +
  ylab("Percent") +
  xlab("Sample")

dir.create("./Figures/Exploration_reviewers/Microbiome/026_Contaminantion_filter", showWarnings = FALSE)
png(paste0("./Figures/Exploration_reviewers/Microbiome/026_Contaminantion_filter/", Tissue, "_Stacked_barchart_contamination_S7_mixed_and_likely.png"),
    res = 600, units = "in", width = 10, height = 4)
plot(plot)
dev.off()

contam_matrix = Genus_full_abundance[which(rownames(Genus_full_abundance) %in% contam),]
contam_matrix$Genus = NULL
contam_matrix = contam_matrix[, which(substring(colnames(contam_matrix), 1, 3) %in% Cohort_246)]
contam_matrix = data.frame(t(contam_matrix))
contam_matrix$Total = rowSums(contam_matrix)

summary(contam_matrix$Total)

write.csv(contam_matrix, file = "./Analysis/Exploration_reviewers/Microbiome/026_Decontamination_prepare_overview/All_S7_contam_in_all_16S_samples.csv")

contam_matrix$Total = NULL
test = data.frame(t(contam_matrix))
test = test ==0
test = test == FALSE
overview = data.frame(rowSums(test))
summary()


