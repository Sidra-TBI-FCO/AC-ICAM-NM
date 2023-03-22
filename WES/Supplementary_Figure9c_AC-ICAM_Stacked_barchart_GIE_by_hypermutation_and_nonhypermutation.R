
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db

setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.bioconductor.packages = c("ggplot2", "ggpubr", "dplyr")                                                                   
ibiopak(required.bioconductor.packages)

# Load data
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
df = Merged_dataset[, c("Patient_ID", "Immunoediting_score", "Immunoedited", "Mutation_cat", "ICR_cluster")]
df$Mutation_cat = factor(df$Mutation_cat, levels = c("nonhypermutated", "hypermutated"))

# Prepare data
DF1 <- df %>%
  group_by(Mutation_cat, Immunoedited) %>%
  dplyr::summarise(count=n()) %>%
  dplyr::mutate(perc=count/sum(count))


plot = ggplot(DF1, aes(x = Mutation_cat, y =perc*100, fill = Immunoedited)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values= c("immunoedited" = "orange", "less immunoedited" = "darkblue"))

dir.create("./Figures/Exploration_reviewers/WES/AC_ICAM_Stacked_barchart_GIE_cat", showWarnings = FALSE)
pdf("./Figures/Exploration_reviewers/WES/AC_ICAM_Stacked_barchart_GIE_cat/AC_ICAM_Stacked_barchart_GIE_cat.pdf",
   width = 3, height = 4)
plot(plot)
dev.off()

chisq.test(table(df$Immunoedited, df$Mutation_cat))
chisq.test(table(df$Mutation_cat, df$Immunoedited))

# Prepare data
DF2 <- df %>%
  group_by(Immunoedited, Mutation_cat) %>%
  dplyr::summarise(count=n()) %>%
  mutate(perc=count/sum(count))


plot = ggplot(DF2, aes(x = Immunoedited, y =perc*100, fill = Mutation_cat)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values= c("nonhypermutated" = "#A7EABD", "hypermutated" = "#EAAED0"))

dir.create("./Figures/Exploration_reviewers/WES/AC_ICAM_Stacked_barchart_GIE_cat", showWarnings = FALSE)
png("./Figures/Exploration_reviewers/WES/AC_ICAM_Stacked_barchart_GIE_cat/AC_ICAM_Stacked_barchart_hypermutation.png",
    res = 600, units = "in", width = 3, height = 4)
plot(plot)
dev.off()

# Prepare data
df$Mutation_cat_ICR = paste(df$Mutation_cat, df$ICR, sep = "_")

df$Mutation_cat_ICR = factor(df$Mutation_cat_ICR, levels = c("hypermutated_ICR High", "nonhypermutated_ICR High",
                                                             "hypermutated_ICR Medium", "nonhypermutated_ICR Medium",
                                                             "hypermutated_ICR Low","nonhypermutated_ICR Low"))

DF3 <- df %>%
  group_by(Mutation_cat_ICR, Immunoedited) %>%
  dplyr::summarise(count=n()) %>%
  mutate(perc=count/sum(count))

plot = ggplot(DF3, aes(x = Mutation_cat_ICR, y =perc*100, fill = Immunoedited)) + geom_bar(stat="identity") +
  labs(x = "", y = "Percentage", fill = "", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 45,
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 19, colour = "black"),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values= c("immunoedited" = "orange", "less immunoedited" = "darkblue"))

png("./Figures/Exploration_reviewers/WES/AC_ICAM_Stacked_barchart_GIE_cat/AC_ICAM_Stacked_barchart_hypermutation_ICR.png",
    res = 600, units = "in", width = 6, height = 7)
plot(plot)
dev.off()

