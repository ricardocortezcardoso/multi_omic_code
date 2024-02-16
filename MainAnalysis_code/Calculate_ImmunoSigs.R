#Script to infer immuno cell signatures in bulk-tumour gene expression data

#Immuno signatures derived from single-cell RNA sequencing data of ccRCC tumours (doi: 10.1016/j.ccell.2022.11.001)
#Input as RNA file (matrix/dataframe of samples as cols and gene expression values as rows) as Rdata file and output directory and dataset name
#This function used data frame with Functional annotation of immuno cells based on a set of specific genes as in Data/Immuno_Signature_Ref.csv
#Scale gene expression values across samples and sum them up by sample.


Calculate_Immuno<-function(RNA_file, dataset_name, output_directory){
  #load RNA file
  RNA_file1=get(load(RNA_file))
  
  #load Immuno cells reference data frame
  Immuno_ref<-as.data.frame(read_csv("Data/Immuno_Signature_Ref.csv"))#provided in Data/
  
  ##create a vector of immune traits
  list_immune=c('NK','B','CD8','CD4','Myeloid','Endothelial','Fibroblasts','Epithelial','RCC')
  
  immune_sig=data.frame(ID=colnames(RNA_file1))
  
  for(i in list_immune){
    print('Selecting annotations > 75% gene overlap')
    
    eval(parse(text=paste0("Immuno_ref %>% filter(Cluster=='",i,"') %>% group_by(annotation) %>% 
    summarise(Freq=length(which(Presence=='Yes'))/length(Presence)) %>%
    filter(Freq>0.75) %>% as.data.frame()->list_annotation")))
    
    for(j in list_annotation$annotation){
      eval(parse(text=paste0("gene_list=Immuno_ref[Immuno_ref$annotation=='",j,"',2]")))
      eval(parse(text=paste0("data=as.data.frame(t(RNA_file1[gene_list,]))")))
      data=data[ , apply(data, 2, function(x) !any(is.na(x)))]
      data %>% mutate(across(.cols =colnames(data), .fns = scale))->data
      eval(parse(text=paste0("temp=data.frame(ID=rownames(data),",j,"=rowSums(data))")))
      immune_sig=left_join(immune_sig,temp,by='ID')
    }
  }
  
  #save file
  write_csv(immune_sig, paste0(output_directory,dataset_name,'.csv'))
  
}

#Example of usage
#Calculate_Immuno(RNA_file ='../../Final_Analysis_RCC/IARC_Expression151RCC_QC.rda',
#                 output_directory = 'multi_omic_code/',dataset_name = 'test'
#                  )

#Expected output is a data frame in CSV  format with immuno signatures as cols and samples as rows