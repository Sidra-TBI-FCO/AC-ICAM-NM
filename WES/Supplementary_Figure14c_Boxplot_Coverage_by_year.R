
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "reshape2", "dplyr", "ggpubr", "stringr", "paletteer"))

# Load data
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")
b1 = read.table("./WES_HPC/AGE_QC/MultiQC_list2001-2005/multiqc_picard_HsMetrics.txt", sep = '\t', header = T)
b2 = read.table("./WES_HPC/AGE_QC/MultiQC_list2006-2010/multiqc_picard_HsMetrics.txt", sep = '\t', header = T)
b3 = read.table("./WES_HPC/AGE_QC/MultiQC_list2011-2015/multiqc_picard_HsMetrics.txt", sep = '\t', header = T)
b_all = rbind(b1, b2)
b_all = rbind(b_all, b3)

# Age
df = data.frame(Patient_ID = clinical_data$Patient_ID, Year_of_collection = clinical_data$year_of_initial_diagnosis, 
                Mean_Target_Coverage = NA, Median_Target_Coverage = NA)

b_all$Sample = gsub("WES_COAD_LUMC_SIDRA_", "", b_all$Sample)
b_all$Sample = str_pad(b_all$Sample, pad = "0", 4)
b_all$Patient_ID = substring(b_all$Sample, 1, 3)

df$Mean_Target_Coverage = b_all$MEAN_TARGET_COVERAGE[match(df$Patient_ID, b_all$Patient_ID)]
df$Median_Target_Coverage = b_all$MEDIAN_TARGET_COVERAGE[match(df$Patient_ID, b_all$Patient_ID)]

palette = paletteer_c("ggthemes::Blue-Green Sequential", 15)

plot = ggplot(df, aes(x = Year_of_collection, y = Mean_Target_Coverage, fill = Year_of_collection)) +
  scale_fill_manual(values = palette) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.8, width = 0.2) +
  theme_bw() +
  xlab("Year of sample collection") +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.position = "none")

pdf("./Figures/Exploration_reviewers/QC/003_Boxplot_Coverage_by_year/003_Supplementary_Figure14c_Boxplot_Mean_Target_Coverage_by_year.pdf",
    width = 7, height = 3)
plot(plot)  
dev.off()



#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(df, varname="Mean_Target_Coverage", 
                    groupnames="Year_of_collection")
# Convert dose to a factor variable
df2$Year_of_collection=as.factor(df2$Year_of_collection)
head(df2)

plot = ggplot(df2, aes(x = Year_of_collection, y = Mean_Target_Coverage, fill = Year_of_collection)) +
  scale_fill_manual(values = palette) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=Mean_Target_Coverage-sd, ymax=Mean_Target_Coverage+sd), width=.2,
                position=position_dodge(.9)) +
  #geom_boxplot(outlier.shape = NA) +
  #geom_jitter(size = 0.8, width = 0.2) +
  theme_bw() +
  xlab("Year of sample collection") +
  ylab("Mean target coverage") +
  ylim(0,max(df$Mean_Target_Coverage)) +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.position = "none")

dir.create("./Figures/Exploration_reviewers/QC/003_Boxplot_Coverage_by_year", showWarnings = FALSE)
pdf("./Figures/Exploration_reviewers/QC/003_Boxplot_Coverage_by_year/003_Boxplot_Mean_Target_Coverage_by_year_barplots.pdf",
    width = 7, height = 3)
plot(plot)  
dev.off()

df$Year_of_collection = factor(df$Year_of_collection)
df$Year_of_collection = as.numeric(df$Year_of_collection)
cor.test(df$Year_of_collection, df$Mean_Target_Coverage, method = "spearman")

plot = ggplot(df, aes(x = Year_of_collection, y = Median_Target_Coverage, fill = Year_of_collection)) +
  scale_fill_manual(values = palette) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.8, width = 0.2) +
  theme_bw() +
  xlab("Year of sample collection") +
  ylab("Median target coverage") +
  ylim(0,max(df$Median_Target_Coverage)) +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.position = "none")

pdf("./Figures/Exploration_reviewers/QC/003_Boxplot_Coverage_by_year/003_Boxplot_Median_Target_Coverage_by_year.pdf",
    width = 7, height = 3)
plot(plot)  
dev.off()

