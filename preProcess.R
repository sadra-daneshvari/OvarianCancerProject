rm(list=ls())
data_folder = "data"
gene_data_folder = "data/Hub genes CSV"
clinical_data_folder = "data/clinical significances"
preProcessed_results_folder = "preProcessed_data"
dir.create(preProcessed_results_folder, showWarnings = FALSE, recursive = TRUE)
all_genes <- c("SYK.csv","CD74.csv","STAT5A.csv","RHOA.csv","CD33.csv",
               "IL1B.csv","WAS.csv","CD163.csv","TLR4.csv","CD4.csv",
               "BST2.csv","CSF1R.csv","TNF.csv","HLA-A.csv","LCP1.csv",
               "B2M.csv","SPI1.csv","HLA-B.csv","ITGAM.csv","NOD2.csv",
               "CEBPB.csv","LYN.csv","CD38.csv","NFKBIA.csv","RAC2.csv",
               "FCER1G.csv","ITGB2.csv","PLEK.csv")
gene_file_dir <- c()
clinical_sigs <- c("Venous invasion.csv", "Primary therapy outcome success.csv",
                   "New neoplasm event type.csv","Neoplasm histologic grade.csv",
                   "Lymphatic invasion.csv")
clinical_sig_dir <- c()
for (gene in all_genes)
{
  gene_file_dir <- append(gene_file_dir, 
                          paste(gene_data_folder,gene,sep="/"))
}
for (clinical_sig in clinical_sigs)
{
  clinical_sig_dir <- append(clinical_sig_dir, 
                            paste(clinical_data_folder,clinical_sig,sep="/"))
}
# build the WHOLE DATA SCV file:
dir_path <- gene_file_dir[1]
info <- read.csv(dir_path)
info_sorted <- info[order(info$sample),]
write.csv(info_sorted,
          paste(preProcessed_results_folder, "WHOLE_DATA.csv",sep="/"))
whole_data_adrs = paste(preProcessed_results_folder,
                        "WHOLE_DATA.csv",
                        sep = "/")
for (dir_path in gene_file_dir[2:28])
{
  info <- read.csv(dir_path)
  info_sorted <- info[order(info$sample),]
  tmp <- read.csv(whole_data_adrs)
  tmp <- cbind(tmp, info_sorted[2])
  write.csv(tmp, whole_data_adrs)
}
for (dir_path in clinical_sig_dir)
{
  info <- read.csv(dir_path)
  info_sorted <- info[order(info$sample),]
  tmp <- read.csv(whole_data_adrs)
  tmp <- cbind(tmp, info_sorted[2])
  write.csv(tmp, whole_data_adrs)
}
whole_data <- read.csv(whole_data_adrs)
cutoffs <- read.csv(paste(data_folder,"CUTOFF.csv", sep="/"), header=FALSE)
cutoff_matrix <- as.matrix(cutoffs[2])
for (i in seq(1:28))
{
  path_gene <- gene_file_dir[i]
  cutoff <- cutoff_matrix[i]
  gene_expression <- read.csv(path_gene)
  gene_exp_mat <- matrix(unlist(gene_expression), ncol=2)
  gene_exp_mat_sorted <- gene_exp_mat[order(gene_exp_mat[,1]),]
  gene_exps <- matrix(gene_exp_mat_sorted[,2], ncol=1)
  num_gene_exps <- matrix(as.numeric(gene_exps))
  colnames(num_gene_exps) <- c(all_genes[i])
  bin_gene_exps <- c()
  for (e in num_gene_exps)
  {
    if (is.na(e))
    {
      bin_gene_exps <- append(bin_gene_exps,NA)
    }
    else
    {
      ifelse (e >= cutoff, 
              bin_gene_exps <- append(bin_gene_exps, "High"),
              bin_gene_exps <- append(bin_gene_exps, "Low"))
    }
  }
  bin_gene_exps_mat <- matrix(bin_gene_exps)
  colnames(bin_gene_exps_mat) <- c(all_genes[i])
  bin_gene_exps_mat
  
  tmp <- read.csv(whole_data_adrs)
  tmp <- cbind(tmp, bin_gene_exps_mat)
  write.csv(tmp, whole_data_adrs)
}

