#Function in R to select the most informative features correlated with MOFA factors

##Function parameters:
##Rdata_path: The path to the Rdata file containing the MOFA object (MOFAobject_ABC) and the dataframe (sampleSheet in CSV format) with sample ID as rows (named as ID) and MOFA factors (Factors1-10) an covariates as columns.
##output_directory: The directory where the output, specifically the LASSO coefficients file, will be saved.
##MOFA_factors: The MOFA factors to be selected for output. This parameter is used to filter the features correlated with the specified MOFA factors.
##Omic: The type of omic data to be used for running the LASSO regression model. It can be either 'RNA' for transcriptome or 'DNAm' for DNA methylome.
##output as dataframe (CSV format)

# Load required R packages
library(tidyverse)#2.0.0
library(tidymodels)#v1.1.1
library(vip)#v0.4.1

generate_MOFAsignatures <- function(Rdata_path, sampleSheet, output_directory, MOFA_factors, Omic) {
  
  # Load Rdata
  load(Rdata_path)
  
  # Load sampleSheet
  sampleSheet1<-as.data.frame(read_csv(sampleSheet))
  
  # Extract relevant information from MOFAobject
  Omic_name<-c(Omic)
  eval(parse(text=paste0("Omic <- as.data.frame(t(as.matrix(MOFAobject_ABC@data$",Omic,"$group1)))")))
  Omic$ID <- rownames(Omic)
  
  # Genomic features correlated with MOFA_factors that were used to train MOFA object
  
  # 1) Get features correlated with MOFA_factors
  temp_data <- sampleSheet1 %>%
    select(ID, MOFA_factors) %>%
    left_join(Omic) %>% rename('Factor'=MOFA_factors) %>% as.data.frame()
  
  results_corr <- data.frame()#col1=sample(ID), col2=(selected MOFA_factor), col3..(other features)
  
  for (i in seq(3, ncol(temp_data))) {
    eval(parse(text = paste0("temp = cor.test(temp_data$Factor, temp_data[, ", i, "])")))
    temp1 <- data.frame(Factor = MOFA_factors, feature = names(temp_data)[i], estimate = temp$estimate, p = temp$p.value)
    results_corr <- rbind(results_corr, temp1)
  }
  
  # Adjust p-values and filter significant results (fdr<0.05)
  results_corr$p.adj <- p.adjust(results_corr$p, method = 'fdr')
  results_corr$Sig <- ifelse(results_corr$p.adj < 0.05, 'Yes', 'No')
  results_corr$Omic_Layer <- Omic_name
  results_corr <- filter(results_corr,Sig == "Yes")
  
  # 2) Run LASSO models
  Train_set <- temp_data %>%
    select(results_corr$feature, Factor, ID) %>%
    mutate(Factor = Factor / sd(Factor)) %>%
    na.omit()  # Remove rows with NAs
  
  ## 2) Recipe - Creating our model
  Train_rec <- recipe(Factor ~ ., data = Train_set) %>%
    update_role(ID, new_role = 'ID') %>%
    step_zv(all_numeric(), -all_outcomes()) %>%
    step_normalize(all_numeric(), -all_outcomes())
  
  # Make sure that the sample variable will not be read as factors by the model
  Train_prep <- Train_rec %>%
    prep(strings_as_factors = FALSE)
  
  # Tune LASSO parameters: picking the penalty value by resampling and tuning
  set.seed(1987)
  Train_boot <- bootstraps(Train_set, times = 1000)
  
  # Tune penalty
  lasso_spec <- linear_reg(penalty = tune(), mixture = 1) %>%
    set_engine("glmnet")
  
  wf <- workflow() %>% add_recipe(Train_rec)
  
  # Make a grid for the amount of penalties (regularization)
  lambda_grid <- grid_regular(penalty(), levels = 100)
  
  # Parallel process for faster computation
  doParallel::registerDoParallel()
  
  set.seed(1987)
  
  lasso_grid <- tune_grid(
    wf %>% add_model(lasso_spec),
    resamples = Train_boot,
    grid = lambda_grid
  )
  
  # Get the best lambda based on rmse
  lowest_rmse<-lasso_grid %>% select_by_one_std_err(metric = 'rmse',maximaze=F)
  lowest_rsq<-lasso_grid %>% select_by_one_std_err(metric = 'rsq',maximaze=F)
  final_lasso <- finalize_workflow(wf %>% add_model(lasso_spec), lowest_rmse)

  # Visualize LASSO coefficients
  results_lasso <- final_lasso %>%
    fit(Train_set) %>%
    pull_workflow_fit() %>%
    vip(lambda = lowest_rmse$penalty)
  
  results_lasso <- as.data.frame(results_lasso$data)
  results_lasso <- mutate(results_lasso,
                       LASSO_coef = ifelse(Sign == 'NEG', Importance * (-1), Importance),
                       Importance = abs(Importance),
                       Variable = fct_reorder(Variable, Importance),
                       rsq=lowest_rsq$mean,
                       rmse=lowest_rmse$mean,
                       Factor= MOFA_factors)
  
  # Save LASSO coefficients
  write_csv(results_lasso, paste0(output_directory, MOFA_factors, '_signature.csv'))
}

#example of usage
#generate_MOFAsignatures(Rdata_path = './Data/Example_MOFA_input.RData',
#                        sampleSheet = './Data/Example_sampleSheet.csv',
#                        MOFA_factors = 'Factor1',
#                        Omic = 'DNAm',
#                        output_directory = './Data/')




