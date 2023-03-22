

# Setup environment
rm(list = ls())

# Set working directory 
load("~/R.Config.Rdata")
master.location = master.location.db
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))
ipak(c("dplyr", "reshape2", "ggplot2", "ggpubr", "dendextend", "corrplot"))

# Set parameters
variable_of_interest = "ESTIMATEScore" #"StromalScore" #"ImmuneScore" #"ESTIMATEScore"
Random = "" # "Random" or ""
size = 200

# Load data
load("./Processed_Data/RNASeq/Trimmed_p/Normalized/JSREP.filtered.dataset.EDAseq.QN.HPC.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/JSREP_clean_ESTIMATE_scores.Rdata")
load("./Analysis/Trimmed_p/ESTIMATE/TCGA_COAD_whitelist_filtered_ESTIMATE_scores.Rdata")

# Data prep
ESTIMATE = ESTIMATE[which(rownames(ESTIMATE) %in% colnames(RNASeq.QN.LOG2)),]
ESTIMATE = data.frame(ESTIMATE)
ESTIMATE_tcga = data.frame(ESTIMATE_tcga)

# Old
ESTIMATE$cohort = "AC-ICAM"
ESTIMATE_tcga$cohort = "TCGA-COAD"
df_plot = rbind(ESTIMATE, ESTIMATE_tcga)

plot = ggplot(df_plot, aes(x = df_plot[,variable_of_interest], color = cohort)) +
  geom_density(aes(fill = cohort), alpha = .4) +
  #xlab(paste0(variable_of_interest)) +
  scale_color_manual(values = c("#ff8c00", "#006DFF")) +
  scale_fill_manual(values = c("#ff8c00", "#006DFF")) +
  theme_bw() +
  xlim(-5000, 5000) +
  xlab(variable_of_interest) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15))

dir.create("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival", showWarnings = FALSE)
pdf(paste0("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Sup_Figure7_New_density_cures_all.pdf"), height = 4, width = 5)
plot(plot)
dev.off() # PLot density curves before sampling


plot2 = ggplot(df_plot, aes(y = df_plot[,variable_of_interest], x = cohort, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.2) +
  #xlab(paste0(variable_of_interest)) +
  scale_fill_manual(values = c("#ff8c00", "#006DFF")) +
  theme_bw() +
  ylim(-5000, 5000) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("AC-ICAM", "TCGA-COAD"))) +
  xlab("") +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x = element_text(color = "black", size = 15),
        axis.title.y = element_text(color = "black", size = 15),
        legend.position = "none")

pdf(paste0("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Sup_Figure7_New_boxplot_all.pdf"), 
    height = 4, width = 3)
plot(plot2)
dev.off() # plot boxplot before sampling

test = t.test(df_plot$ESTIMATEScore[which(df_plot$cohort == "AC-ICAM")], df_plot$ESTIMATEScore[which(df_plot$cohort == "TCGA-COAD")])
test$p.value

# Sampling from Sidra-LUMC cohort
ESTIMATE_full= ESTIMATE

i=81
for(i in 1:100){
  set.seed(i)
  
  # Obtain the density of the TCGA at the length of the model.
  dens.tcga <- density(ESTIMATE_tcga$ESTIMATEScore, n=size)
  
  predict_density = approxfun(dens.tcga) #function that approximates dens.tcga
  #sample points from Model with probability distr. of dens.obs
  Sample_from_SILU <- sample(rownames(ESTIMATE_full), size = size, replace=FALSE, prob=predict_density(ESTIMATE_full$ESTIMATEScore))
  Random_Sample_from_SILU <- sample(rownames(ESTIMATE_full), size = size, replace=FALSE)
  
  if(Random == "Random"){
    ESTIMATE = ESTIMATE_full[Random_Sample_from_SILU,]
  }else{
    ESTIMATE = ESTIMATE_full[Sample_from_SILU,]
  }
  
  
  # New
  ESTIMATE$cohort = "AC-ICAM"
  ESTIMATE_tcga$cohort = "TCGA-COAD"
  df_plot = rbind(ESTIMATE, ESTIMATE_tcga)
  
  plot = ggplot(df_plot, aes(x = df_plot[,variable_of_interest], color = cohort)) +
    geom_density(aes(fill = cohort), alpha = .4) +
    #xlab(paste0(variable_of_interest)) +
    scale_color_manual(values = c("#ff8c00", "#006DFF")) +
    scale_fill_manual(values = c("#ff8c00", "#006DFF")) +
    theme_bw() +
    xlab(variable_of_interest) +
    xlim(-5000, 5000) +
    theme(axis.text.x = element_text(color = "black", size = 11),
          axis.text.y = element_text(color = "black", size = 11),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15))
  
  
  dir.create("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival", showWarnings = FALSE)
  pdf(paste0("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Sup_Figure7_New_density_curves_", Random, "_set_seed_", i,".pdf"), height = 4, width = 5)
  plot(plot)
  dev.off()
  
  plot2 = ggplot(df_plot, aes(y = df_plot[,variable_of_interest], x = cohort, fill = cohort)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.2) +
    #xlab(paste0(variable_of_interest)) +
    scale_fill_manual(values = c("#ff8c00", "#006DFF")) +
    theme_bw() +
    stat_compare_means(method = "t.test", label = "p.signif",
                       comparisons = list(c("AC-ICAM", "TCGA-COAD"))) +
    xlab(paste0(variable_of_interest)) +
    ylim(-5000, 5000) +
    theme(axis.text.x = element_text(color = "black", size = 11),
          axis.text.y = element_text(color = "black", size = 11),
          axis.title.x = element_text(color = "black", size = 15),
          axis.title.y = element_text(color = "black", size = 15),
          legend.position = "none")
  
  pdf(paste0("./Figures/Exploration_reviewers/Purity_cutoff_ICR_Survival/Sup_Figure7_New_boxplot_", Random, "_set_seed_", i,".pdf"), 
      height = 4, width = 3)
  plot(plot2)
  dev.off()
  
  dir.create("./Analysis/Exploration_reviewers/Purity_cutoff_ICR_Survival", showWarnings = FALSE)
  save(ESTIMATE, file = paste0("./Analysis/Exploration_reviewers/Purity_cutoff_ICR_Survival/New_ESTIMATE_", Random,
                               "set_seed_", i, ".Rdata"))
  
}

#set.seed(12) #12




