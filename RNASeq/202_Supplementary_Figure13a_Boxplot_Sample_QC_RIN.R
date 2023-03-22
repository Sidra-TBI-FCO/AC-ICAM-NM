
# Set-up environment
rm(list = ls())

load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

ipak(c("ggplot2", "reshape2", "dplyr", "ggpubr", "stringr", "paletteer"))

# Load data
RIN_RNASeq = read.csv("./Sample_QC/Sheet2_JSREP RNA TUMOUR QC REPORT Sept 2022.csv", stringsAsFactors = FALSE)
load("./Processed_Data/Survival Data/clinical_data_348_patients.Rdata")

# Age
df = data.frame(Patient_ID = clinical_data$Patient_ID, Year_of_collection = clinical_data$year_of_initial_diagnosis, RIN = NA)

RIN_RNASeq$Sample.Name = str_pad(RIN_RNASeq$Sample.Name, pad = "0", 4)
RIN_RNASeq = RIN_RNASeq[-(grep("LM", RIN_RNASeq$Sample.Name)),]
RIN_RNASeq = RIN_RNASeq[-(grep("B", RIN_RNASeq$Sample.Name)),]
RIN_RNASeq = RIN_RNASeq[-(grep("C", RIN_RNASeq$Sample.Name)),]
df$RIN = RIN_RNASeq$RIN[match(df$Patient_ID, substring(RIN_RNASeq$Sample.Name, 1, 3))]
df$RIN = as.numeric(df$RIN)

palette = paletteer_c("ggthemes::Blue-Green Sequential", 15)

# RIN by diagnosis year
plot = ggplot(df, aes(x = Year_of_collection, y = RIN, fill = Year_of_collection)) +
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

pdf("./Figures/Exploration_reviewers/QC/002_RIN_by_year_of_collection/002_Supplementary_Figure13a_RIN_by_year_of_collection_barplot.pdf",
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

df2 <- data_summary(df, varname="RIN", 
                    groupnames="Year_of_collection")
# Convert dose to a factor variable
df2$Year_of_collection=as.factor(df2$Year_of_collection)
head(df2)

# RIN by diagnosis year
plot = ggplot(df2, aes(x = Year_of_collection, y = RIN, fill = Year_of_collection)) +
  scale_fill_manual(values = palette) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=RIN-sd, ymax=RIN+sd), width=.2,
                position=position_dodge(.9)) +
  #geom_jitter(size = 0.8, width = 0.2) +
  theme_bw() +
  xlab("Year of sample collection") +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.position = "none")

dir.create("./Figures/Exploration_reviewers/QC/002_RIN_by_year_of_collection", showWarnings = FALSE)
pdf("./Figures/Exploration_reviewers/QC/002_RIN_by_year_of_collection/002_RIN_by_year_of_collection_barplot.pdf",
    width = 7, height = 3)
plot(plot)  
dev.off()



