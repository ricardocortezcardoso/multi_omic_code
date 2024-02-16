####################################################################################################
##########################Pathway analysis of Genes correlated with latent factors################## 
#####################################################################################################
#devtools::install_github("GSEA-MSigDB/GSEA_R")
#devtools::install_github("ctlab/fgsea")
#GSEA code 
#GO_file is annotation file; 
#named numeric vector containing gene symbols as name, and correlation/estimate as numbers
#this function performs a pathway analysis based on a set of genes correlated with outcome of interest (e.g., MOFA factors)
#Inputs are 1) correlation file as in Data/Example_Correlation_Factor3_RNA.csv; 2)full path to output directoty 3) outcome name (e.g., Factor3)

set.seed(54321)
library(tidyverse)#v2.0.0
library(fgsea)#v1.27.1
library(dplyr)#v1.1.3

pathway_analysis<-function(Correlation_file, output_folder, Outcome){
  
  #1) Upload GO terms to Hallmarks curate from Chen et al., 2021 (https://doi.org/10.1186/s12859-021-04105-8)
  GO_Hallmark<-read_csv("Data/Hallmarks_Cancer_GO.csv")
  
  #2) GO terms to Pathways
  GO_Path<-read_csv("Data/GOTerm_dictionary.csv")
  
  #3) Upload files with correlations between Factor and Gene expression as in ./Data/Example_Correlation_Factor3_RNA.csv
  results_corr<-as.data.frame(read_csv(Correlation_file))
  
  #4) Select the top 500 associations that passed p.adj (FDR) <0.05
  results_corr %>% filter(Sig=='Yes') %>% arrange(p) %>% dplyr::slice(1:500)->results_corr
  
  geneList = as.numeric(results_corr[,3])#third col contains the rho estimates (Perason's correlations between Gene and Factor)
  names(geneList)= as.character(results_corr[,2]) # Vector with Gene symbols
  
  #5) Run fgsea
  #5.1) Open annotation file (https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/msigdb.v2022.1.Hs.symbols.gmt)
  
  myGO = fgsea::gmtPathways("Data/msigdb.v2022.1.Hs.symbols.gmt")
  
  #5.2) Run fgsea 
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = geneList,
                        minSize=15, ## minimum gene set size
                        maxSize=5000,
                        eps=0.0) %>% ## maximum gene set size
    as.data.frame() %>% 
    dplyr::filter(padj < 0.05) %>% 
    arrange(desc(NES))
  
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = geneList)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  fgRes$Factor=Outcome
  fgRes=left_join(fgRes,GO_Path,by='pathway') #get pathway enrichment (standard)
  fgRes=left_join(fgRes,GO_Hallmark,by='GO') # get Hallmarks enrichment
  
  #save results
  write_csv(fgRes,paste0(output_folder, Outcome,'.csv'))
  
}

#Example of usage
#pathway_analysis(Correlation_file = './Data/Example_Correlation_Factor3_RNA.csv',
#                 output_folder = './Data/',
#                 Outcome = 'Factor3'
#)


#Expected output: a data frame in CSV format with pathway, pvalues, Enrichment scores, Direction of association, Outcome, GO, Hallmark of cancer