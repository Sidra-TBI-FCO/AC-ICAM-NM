
# Set-up environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

# load library
required.bioconductor.packages = c("ggplot2", "ggpubr")                                                                   
ibiopak(required.bioconductor.packages)

#load("./Analysis/TCR/0_1_Annotated_data/0_1_Annotated_TCR_Overview.Rdata")
load("./Analysis/WES/022b_Expected_versus_Observed_Neoantigens/022c_Merged_dataset.Rdata")
load("./Analysis/WES/022f_IES_categories/022f_IES_df.Rdata")
MiXCR_orig = read.csv("./Processed_Data/MiXCR/RNASeq_Mixcr_results_19_September_2020/Diversity.csv", stringsAsFactors = FALSE)

#colnames(Merged_dataset)

Merged_dataset$IES = IES_df$IES[match(Merged_dataset$Patient_ID, IES_df$Patient_ID)]
colnames(Merged_dataset)

# Prepare MiXCR data
MiXCR = MiXCR_orig
colnames(MiXCR)[1] = "Sample_ID"
MiXCR$Patient_ID = substring(MiXCR$Sample_ID, 1, 3)
MiXCR$Tissue = substring(MiXCR$Sample_ID, 4, 6)

MiXCR = MiXCR[which(MiXCR$Tissue == "T-P"),] # Filter only tumor samples
MiXCR = MiXCR[which(substring(MiXCR$Sample_ID, 7, 10) == "_TRB"),]
MiXCR = MiXCR[-which(is.na(MiXCR$Pielou)),]
MiXCR$TCR_Clonality = 1 - MiXCR$Pielou

MiXCR = MiXCR[-which(MiXCR$TCR_Clonality == -Inf),]

# Overview subgroups
Merged_dataset$TCR_MiXCR_clonality = MiXCR$TCR_Clonality[match(Merged_dataset$Patient_ID,
                                                               MiXCR$Patient_ID)]
table(Merged_dataset$Immunoedited)

# Darawans plot 

df.plot = Merged_dataset[,c(1,65,66)]

df.plot = df.plot[!is.na(df.plot$IES),]
df.plot = df.plot[!is.na(df.plot$TCR_MiXCR_clonality),]

df.plot$IES.new = df.plot$IES

df.plot$IES = as.numeric(factor(df.plot$IES, levels = c("IES1", "IES2", "IES3", "IES4")))

df.plot$TCR_productive_clonality = Merged_dataset$TCR_productive_clonality[match(df.plot$Patient_ID, Merged_dataset$Patient_ID)]

#dir.create("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal")
svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/geom_smooth_IES_TCR_MiXCR_included.svg"), width = 3, height = 3)

plot.df = ggplot(df.plot, aes(x = IES, y = TCR_MiXCR_clonality)) +
  geom_smooth(method=loess, se = TRUE) +
  #scale_color_manual(values = c("orange", "yellow", "purple", "brown"))+
  #scale_fill_manual(values = c("red", "blue", "green", "black"))+
  #scale_x_discrete(name = "IES", limits=c("IES1","IES2",
  #"IES3", "IES4")) +
  scale_y_continuous(name = "TCR_MiXCR_clonality") +
  xlab("IES") +
  guides(color=guide_legend(ncol=1))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=15,hjust=0.5),
        legend.text = element_text(size=15),
        axis.text.x = element_text(colour="black",size=15,angle=360,hjust=1,vjust=1,face="plain"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=15,face="plain"),
        #axis.title.x = element_text(colour="black",size=30,face="plain"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black",size=15,face="plain"),
        #legend.title = element_blank(),
        legend.position = "right")

print(plot.df)
dev.off()

svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/geom_smooth_IES_TCR_productive_clonality.svg"), width = 3, height = 3)

plot.df = ggplot(df.plot, aes(x = IES, y = TCR_productive_clonality)) +
  geom_smooth(method=loess, se = TRUE) +
  #scale_color_manual(values = c("orange", "yellow", "purple", "brown"))+
  #scale_fill_manual(values = c("red", "blue", "green", "black"))+
  #scale_x_discrete(name = "IES", limits=c("IES1","IES2",
  #"IES3", "IES4")) +
  scale_y_continuous(name = "TCR_productive_clonality") +
  xlab("IES") +
  guides(color=guide_legend(ncol=1))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=15,hjust=0.5),
        legend.text = element_text(size=15),
        axis.text.x = element_text(colour="black",size=15,angle=360,hjust=1,vjust=1,face="plain"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=15,face="plain"),
        #axis.title.x = element_text(colour="black",size=30,face="plain"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black",size=15,face="plain"),
        #legend.title = element_blank(),
        legend.position = "right")

print(plot.df)
dev.off()

pdf(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/geom_smooth_IES_TCR_productive_clonality_v1.pdf"), width = 3, height = 3)

plot.df = ggplot(df.plot, aes(x = IES, y = TCR_productive_clonality)) +
  geom_smooth(method=loess, se = TRUE) +
  #scale_color_manual(values = c("orange", "yellow", "purple", "brown"))+
  #scale_fill_manual(values = c("red", "blue", "green", "black"))+
  #scale_x_discrete(name = "IES", limits=c("IES1","IES2",
  #"IES3", "IES4")) +
  scale_y_continuous(name = "TCR_productive_clonality") +
  xlab("IES") +
  guides(color=guide_legend(ncol=1))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", size=15,hjust=0.5),
        legend.text = element_text(size=15),
        axis.text.x = element_text(colour="black",size=15,angle=360,hjust=1,vjust=1,face="plain"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=15,face="plain"),
        #axis.title.x = element_text(colour="black",size=30,face="plain"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black",size=15,face="plain"),
        #legend.title = element_blank(),
        legend.position = "right")

print(plot.df)
dev.off()

####################################################################################
# 
# Merged_dataset$IES.new = NA
# 
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR High" & Merged_dataset$Immunoedited == "immunoedited")] = "IES6"
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR High" & Merged_dataset$Immunoedited == "less immunoedited")] = "IES5"
# 
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR Medium" & Merged_dataset$Immunoedited == "immunoedited")] = "IES4"
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR Medium" & Merged_dataset$Immunoedited == "less immunoedited")] = "IES3"
# 
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR Low" & Merged_dataset$Immunoedited == "immunoedited")] = "IES2"
# Merged_dataset$IES.new[which(Merged_dataset$ICR_cluster == "ICR Low" & Merged_dataset$Immunoedited == "less immunoedited")] = "IES1"
# 
# df.plot = Merged_dataset[,c(1,66,67)]
# 
# df.plot = df.plot[!is.na(df.plot$IES.new),]
# df.plot = df.plot[!is.na(df.plot$TCR_MiXCR_clonality),]
# 
# df.plot$IES = df.plot$IES.new
# 
# df.plot$IES = as.numeric(factor(df.plot$IES, levels = c("IES1", "IES2", "IES3", "IES4", "IES5", "IES6")))
# 
# svg(paste0("./Figures/Exploration_reviewers/Microbiome/tumor_only_246/rebuttal/geom_smooth_IES_TCR_MiXCR_included_medium.svg"), width = 3, height = 3)
# 
# plot.df = ggplot(df.plot, aes(x = IES, y = TCR_MiXCR_clonality)) +
#   geom_smooth(method=loess, se = TRUE) +
#   #scale_color_manual(values = c("orange", "yellow", "purple", "brown"))+
#   #scale_fill_manual(values = c("red", "blue", "green", "black"))+
#   #scale_x_discrete(name = "IES", limits=c("IES1","IES2",
#   #"IES3", "IES4")) +
#   scale_y_continuous(name = "TCR_MiXCR_clonality") +
#   xlab("IES") +
#   guides(color=guide_legend(ncol=1))+
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(face="bold", size=15,hjust=0.5),
#         legend.text = element_text(size=15),
#         axis.text.x = element_text(colour="black",size=15,angle=360,hjust=1,vjust=1,face="plain"),
#         #axis.text.x = element_blank(),
#         axis.text.y = element_text(colour="black",size=15,face="plain"),
#         #axis.title.x = element_text(colour="black",size=30,face="plain"),
#         #axis.title.x = element_blank(),
#         axis.title.y = element_text(colour="black",size=15,face="plain"),
#         #legend.title = element_blank(),
#         legend.position = "right")
# 
# print(plot.df)
# dev.off()
