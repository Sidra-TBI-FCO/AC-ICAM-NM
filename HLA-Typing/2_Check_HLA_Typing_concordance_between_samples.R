

## Check matching of HLA type by sample

# Set-up environment
rm(list = ls())
load("~/R.Config.Rdata")
setwd(paste0(master.location, "/TBI-LAB - Project - COAD SIDRA-LUMC/NGS_Data"))
source(paste0(toolbox.path, "/R scripts/ipak.function.R"))

required.packages = c("stringr")
ipak(required.packages)

# Load data
load("./Analysis/HLA Typing/Optitype_merged_result/Results_Optitype_RNASeq_WES_by_Patient_ID.Rdata")
all_results_df$Patient_ID <- as.character(all_results_df$Patient_ID)

# Add column to check if the columns match
all_results_df$concordance = NA

for(i in 1:nrow(all_results_df)){
  test.data = all_results_df[i,c(4:8)]            #select data
  test.data = test.data[-which(is.na(test.data))] #remove NA data
  #print(test.data)
  all_results_df$concordance[i] = all(sapply(as.list(test.data),function(x) x==test.data[1]))
}

dir.create("./Analysis/HLA Typing/Optitype_concordance", showWarnings = FALSE)
save(all_results_df ,file = "./Analysis/HLA Typing/Optitype_concordance/Results_Optitype_concordance.Rdata")
write.csv(all_results_df, file = "./Analysis/HLA Typing/Optitype_concordance/Results_Optitype_concordance.csv")

# check diconcordant patients 
discon.patients <- unique(all_results_df$Patient_ID[which(all_results_df$concordance == FALSE)]) #46
discon_results_df <- all_results_df[all_results_df$Patient_ID %in% discon.patients,]
discon_details_df <- data.frame(patient_ID = discon.patients,
          A.match.T_RNA_T_WES = NA,
          B.match.T_RNA_T_WES = NA,
          C.match.T_RNA_T_WES = NA,
          A.match.T_WES_N_WES = NA,
          B.match.T_WES_N_WES = NA,
          C.match.T_WES_N_WES = NA,
          A.match.T_RNA_N_WES = NA,
          B.match.T_RNA_N_WES = NA,
          C.match.T_RNA_N_WES = NA)
i="110"
for(i in discon.patients){
  HLA.data <- discon_results_df[discon_results_df$Patient_ID==i,]
  for (j in c("A","B","C")) {
    assign(paste0(j,".match.T_RNA_T_WES"),HLA.data[HLA.data$Locus==j,"Alleles_T_RNASeq"] == HLA.data[HLA.data$Locus==j,"Alleles_T_WES"])
    assign(paste0(j,".match.T_WES_N_WES"),HLA.data[HLA.data$Locus==j,"Alleles_T_WES"] == HLA.data[HLA.data$Locus==j,"Alleles_N_WES"])
    assign(paste0(j,".match.T_RNA_N_WES"),HLA.data[HLA.data$Locus==j,"Alleles_T_RNASeq"] == HLA.data[HLA.data$Locus==j,"Alleles_N_WES"])
  }
  discon_details_df[discon_details_df$patient_ID == i,c(2:ncol(discon_details_df))] <- c(A.match.T_RNA_T_WES,B.match.T_RNA_T_WES,C.match.T_RNA_T_WES,
                                                                                         A.match.T_WES_N_WES,B.match.T_WES_N_WES,C.match.T_WES_N_WES,
                                                                                         A.match.T_RNA_N_WES,B.match.T_RNA_N_WES,C.match.T_RNA_N_WES)
}
discon_details_df[discon_details_df$patient_ID == "110",c(2:ncol(discon_details_df))] <- c(rep(FALSE,9))

#conclusion
discon_details_df$TRTW.mismatch = rowSums(!discon_details_df[,c("A.match.T_RNA_T_WES","B.match.T_RNA_T_WES","C.match.T_RNA_T_WES")])
discon_details_df$TWNW.mismatch = rowSums(!discon_details_df[,c("A.match.T_WES_N_WES","B.match.T_WES_N_WES","C.match.T_WES_N_WES")])
discon_details_df$TRNW.mismatch = rowSums(!discon_details_df[,c("A.match.T_RNA_N_WES","B.match.T_RNA_N_WES","C.match.T_RNA_N_WES")])
discon_details_df$conclusion = NA

discon_details_df$conclusion[discon_details_df$TWNW.mismatch > 0 & discon_details_df$TRNW.mismatch > 0] <- "double mismatch"
discon_details_df$conclusion[discon_details_df$TWNW.mismatch == 0 & discon_details_df$TRNW.mismatch > 0] <- "TRNW mismatch"
discon_details_df$conclusion[discon_details_df$TWNW.mismatch > 0 & discon_details_df$TRNW.mismatch == 0] <- "TWNW mismatch"

discon_details_df$conclusion[discon_details_df$patient_ID == "110"] <- "prim-meta mismatch"

write.csv(discon_details_df,file="./Analysis/HLA Typing/Optitype_concordance/disconcordancy_details.csv")  

