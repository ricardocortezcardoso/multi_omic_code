##Inferring MOFA factors in external datasets

library(tidyverse)
library(methylkey)
library(SummarizedExperiment)

infer_MOFA_factors <- function(path_to_LASSO_coefficients, path_to_Omic_file, dataset_name, output_directory, MOFA_factor) {
  # Open LASSO coefficients for each MOFA factor
  LASSO_coef <- as.data.frame(read_csv(path_to_LASSO_coefficients))
  
  # Get DNA methylome (m-values) or RNA (log2-voom transformed if RNA-seq) in TXT format, rows as feature, cols as samples
  Omic_file <- bigreadr::fread2(path_to_Omic_file)
  
  # Inferring signatures
  signature_file <- data.frame(ID = colnames(Omic_file))
  
  # Checking how many LASSO features overlap with Omic_file
  temp_data <- as.data.frame(t(Omic_file[rownames(Omic_file) %in% LASSO_coef$Variable, ]))
  if (nrow(LASSO_coef) == ncol(temp_data)) {
    print("All features found")
  }
  if (ncol(temp_data) / nrow(LASSO_coef) < 0.75) {
    print(paste0('Less than 75% found (N=', ncol(temp_data), ')'))
  }
  LASSO_coef <- filter(LASSO_coef, Variable %in% colnames(temp_data))
  
  # Prepare file
  temp_data %>% mutate(across(.cols = names(temp_data)[1:ncol(temp_data)], .fns = scale)) %>% as.matrix() -> temp_data
  temp_data[!is.finite(temp_data)] <- 0
  temp_data <- as.data.frame(temp_data)
  temp_data$ID <- rownames(temp_data)
  
  # Multiply values by respective LASSO coefs
  for (p in LASSO_coef$Variable) {
    eval(parse(text = paste0("temp_data <- mutate(temp_data, ", p, "=", p, "*LASSO_coef[which(LASSO_coef$Variable=='", p, "'),4])")))
  }
  
  # Sum up values across features by sample
  data_matrix <- t(temp_data[, 1:10])
  eval(parse(text = paste0("temp = data.frame(ID = colnames(data_matrix), Factor", j, " = colSums(data_matrix))")))
  signature_file <- left_join(signature_file, temp)
  
  # calculate variance across samples in Omic file explained by MOFA factor
  Omic_file %>% dplyr::select(signature_file$ID)->Omic_file
  all(colnames(Omic_file)==signature_file$ID)#TRUE
  rownames(signature_file)=signature_file$ID
  R2_met.pro = colMeans(t(cor(signature_file[,which(colnames(signature_file)==MOFA_factor)],t(Omic_file),use="pairwise.complete.obs") )**2,na.rm=T)  
  result_file=rbind(result_file,data.frame(r2=R2_met.pro,cohort=dataset_name))
  
  #save
  write_csv(signature_file, paste0(output_directory, dataset_name, '_', MOFA_factor, '.csv'))
  write_csv(result_file, paste0(output_directory, 'R2_',dataset_name, '_', MOFA_factor, '.csv'))
}

# Example usage:
# infer_MOFA_factors("path/to/LASSO_coefficients.csv", "path/to/Omic_file.txt", "dataset_name" ,"output/directory/", "MOFA_factor_name")


