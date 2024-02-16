##R script to calculate Epigenetic clocks 
##dnaMethyAge from Github page (https://github.com/yiluyucheng/dnaMethyAge)
##devtools::install_github("yiluyucheng/dnaMethyAge")
#This function uses dnaMethyAge R package to infer epigenetic age in a dataset. check the availableClock() to see which clocks are available 
#This function uses normalised (Funnrom) beta values (output from pre_processing_DNAMethyl.Rmd) and sample_file (example in Data/Example_sampleSheet.csv), and clock name (vector with clock names as strings)

library(dnaMethyAge)#v0.1.0

availableClock()#check the available clocks 

Epigenetic_clock<-function(DNA_file,sampleSheet, clock_name, output_dir, dataset_name){
  
  #load DNAm_file
  load(DNA_file)
  
  #load sampleSheet file
  sampleSheet1<-as.data.frame(read_csv(sampleSheet))
  sampleSheet1<-filter(sampleSheet1, Basename %in% colnames(betas))#assuming same IDs
  
  #create a dataframe with Sample, Age and Sex cols (mandatory colnames for the dnaMethyAge)
  info=data.frame(Sample=sampleSheet1$Basename,Age=sampleSheet1$age_diag, Sex=sampleSheet1$sex)
  
  #calculate epigenetic age by clock
  
  clock_res=data.frame(Sample=colnames(betas))
  
  for (i in clock_name){
    temp <- dnaMethyAge::methyAge(betas, clock=i, age_info=info, fit_method='Linear', do_plot=F)
    temp %>% dplyr::select(Sample,mAge,Age_Acceleration)->temp
    names(temp)[2:3]=c(i,paste0('adj_',i))
    clock_res=left_join(clock_res,temp,by='Sample')
  }
  
  write_csv(clock_res,paste0(output_dir,"EpigeneticClocks_",dataset_name,'.csv'))
    
}

#Example of usage
#Epigenetic_clock(DNA_file ="/data/geniluc/work/cortezr/RCC_Mutographs/RCC_121_Funnorm_2023-09-04.rda",
#                 output_dir ="multi_omic_code/",sampleSheet = "multi_omic_code/Data/Example_sampleSheet.csv",
#                 clock_name = c("DunedinPACE","epiTOC2"),
#                 dataset_name = "test"
#)

#Expected output: a data frame in CSV format with epigenetic clocks (unadjusted and age-adjusted) as cols and samples as rows






