rm(list=ls())
preProcessed_results_folder = "preProcessed_data"
whole_data_adrs = paste(preProcessed_results_folder,
                        "WHOLE_DATA.csv",
                        sep="/")
results_adrs = "results"
dir.create(results_adrs, showWarnings = FALSE, recursive = TRUE)
final_data <- read.csv(whole_data_adrs)
final_data <- na.omit(final_data)
final_data = final_data[,c("sample",
                           "venous_invasion",
                           "primary_therapy_outcome_success",
                           "new_neoplasm_event_type",
                           "neoplasm_histologic_grade",
                           "lymphatic_invasion",
                           "SYK.csv",
                           "CD74.csv",
                           "STAT5A.csv",
                           "RHOA.csv",
                           "CD33.csv",
                           "IL1B.csv",
                           "WAS.csv",
                           "CD163.csv",
                           "TLR4.csv",
                           "CD4.csv",
                           "BST2.csv",
                           "CSF1R.csv",
                           "TNF.csv",
                           "HLA-A.csv",
                           "LCP1.csv",
                           "B2M.csv",
                           "SPI1.csv",
                           "HLA-B.csv",
                           "ITGAM.csv",
                           "NOD2.csv",
                           "CEBPB.csv",
                           "LYN.csv",
                           "CD38.csv",
                           "NFKBIA.csv",
                           "RAC2.csv",
                           "FCER1G.csv",
                           "ITGB2.csv",
                           "PLEK.csv")]
# data of 308 patients without NA values built!
# NOW altering neoplasm_histologic_grade; dividing into 2 groups: g1+g2 and g3+g4 (ignore GX, GB)
altered_neoplasm_grade <- c()
for (i in final_data$neoplasm_histologic_grade)
{
  if (i == "G3" || i == "G4")
  {
    altered_neoplasm_grade<-append(altered_neoplasm_grade, "G3+G4")
  }
  else
  {
    ifelse(i == "G2" || i == "G1", altered_neoplasm_grade<-append(altered_neoplasm_grade, "G1+G2"),
           altered_neoplasm_grade<-append(altered_neoplasm_grade, NA))
  }
}
altered_neoplasm_grade <- matrix(altered_neoplasm_grade)
final_data$neoplasm_histologic_grade <- altered_neoplasm_grade
write.csv(final_data,
          paste(results_adrs, "final_data.csv", sep="/"))
# finally built the final data for #308 patients without NA values in their gene expressions
# NOW: chi-sq test on all genes for all attributes and building p-value csv file
chi_matrix <- matrix(NA,nrow = 5, ncol = 28)
colnames(chi_matrix) <- c("SYK.csv",
                          "CD74.csv",
                          "STAT5A.csv",
                          "RHOA.csv",
                          "CD33.csv",
                          "IL1B.csv",
                          "WAS.csv",
                          "CD163.csv",
                          "TLR4.csv",
                          "CD4.csv",
                          "BST2.csv",
                          "CSF1R.csv",
                          "TNF.csv",
                          "HLA.A.csv",
                          "LCP1.csv",
                          "B2M.csv",
                          "SPI1.csv",
                          "HLA.B.csv",
                          "ITGAM.csv",
                          "NOD2.csv",
                          "CEBPB.csv",
                          "LYN.csv",
                          "CD38.csv",
                          "NFKBIA.csv",
                          "RAC2.csv",
                          "FCER1G.csv",
                          "ITGB2.csv",
                          "PLEK.csv")
rownames(chi_matrix) <- c("venous_invasion",
                          "primary_therapy_outcome_success",
                          "new_neoplasm_event_type",
                          "neoplasm_histologic_grade",
                          "lymphatic_invasion")
for (gene in colnames(chi_matrix))
{
  for (factor in rownames(chi_matrix))
  {
    if(length(unique(final_data[,gene])) > 1)
    {
      chi_test = chisq.test(final_data[, factor], final_data[, gene])
      chi_matrix[factor, gene] <- chi_test$p.value
    }
  }
}
write.csv(chi_matrix,
          paste(results_adrs,"p_values.csv",sep="/"))
p_adj_mat = matrix(p.adjust(chi_matrix), nrow(chi_matrix), ncol(chi_matrix))
rownames(p_adj_mat) <- rownames(chi_matrix)
colnames(p_adj_mat) <- colnames(chi_matrix)
write.csv(p_adj_mat, 
          paste(results_adrs, "pVal_adj.csv", sep="/"))
  
