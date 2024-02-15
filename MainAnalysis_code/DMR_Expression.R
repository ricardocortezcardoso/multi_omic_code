#Identifying DMRs using DNA methylome and how they correlate to gene expression levels
#BiocManager::install("DMRcate")
#devtools::install_github("perishky/dmrff")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

##DNAm_file: path to rda file with processed (betas and mval) and and RGset methylation files, output from pre_processing_DNAMethyl.Rmd
##RNA_file: path to rda file with processed gene expression data (genes as rows, sample as cols)
##sampleSheet: info file wit patients data (features as cols, samples as rows), including sequencing (Basename) and epi names (ID) as CSV format


#load packages
library(tidyverse)#v2.0.0
library(data.table)#v1.14.8
library(dmrff)#v1.1.1
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)#v0.6.0
library(limma)#v3.50.3


calculate_DMR <- function(DNAm_file, RNA_file, sampleSheet, MOFA_factors, output_directory) {
  
  print("loading sampleSheet file")
  
  sampleSheet1<-as.data.frame(read_csv(sampleSheet))
  sampleSheet1<-filter(sampleSheet1,!is.na(Basename))#exclude samples without DNA methylation
  
  print("load normalized DNAm files")
  
  load(DNAm_file)
  if(all(colnames(betas) %in% sampleSheet1$Basename)){
    betas=as.data.frame(betas)
    data.table::setnames(betas, old=sampleSheet1$Basename, new=sampleSheet1$ID,skip_absent=TRUE)
    betas=as.matrix(betas)
    
    mval=as.data.frame(mval)
    data.table::setnames(mval, old=sampleSheet1$Basename, new=sampleSheet1$ID,skip_absent=TRUE)
    mval=as.matrix(mval)
  }
  
  print("get annotation file")
  
  Full_ann=as.data.frame(minfi::getAnnotation(RGset, lociNames = NULL,
                                              orderByLocation = FALSE, dropNonMapping = FALSE))
  Full_ann %>%  
    mutate(UCSC_RefGene_Name=gsub(";.*",'',UCSC_RefGene_Name),
           UCSC_RefGene_Group=gsub(";.*",'',UCSC_RefGene_Group),
           enhancer=ifelse(nchar(Phantom4_Enhancers)>1 | nchar(Phantom5_Enhancers)>1 | nchar(X450k_Enhancer)>1 ,'Yes','No'),
           OpenChromatin=ifelse(!OpenChromatin_Evidence_Count=='','Open_Chromatin',''),
           TFBS=ifelse(!TFBS_Evidence_Count=='','TFBS',''))->Full_ann
  
  print("load RNA file")
  RNA_file1=get(load(RNA_file))
  RNA_file1=RNA_file1[,colnames(RNA_file1) %in% colnames(mval)]
  
  #final QC for methylation data
  noSNPs <- DMRcate::rmSNPandCH(betas, dist=2, mafcut=0.05)
  nrow(noSNPs)
  mval <- mval[rownames(noSNPs),] 
  colnames(noSNPs) <- colnames(mval)
  
  print("running limma")
  
  sampleSheet1<-dplyr::rename(sampleSheet1,'Factor'=MOFA_factors)
  
  design <- model.matrix(~ Factor+age_diag+sex,data = sampleSheet1)#covariates can be added here. 
  fit <- limma::lmFit(mval, design,method = 'ls')
  fit <- limma::eBayes(fit)
  stats<-data.frame(estimate=fit$coefficients[,"Factor"],
                    se=sqrt(fit$s2.post) * fit$stdev.unscaled[,"Factor"],
                    p.value=fit$p.value[,"Factor"])
  
  #creating annotation file
  common <- intersect(rownames(mval), rownames(Full_ann))
  Full_ann <- Full_ann[match(common, rownames(Full_ann)),]
  stats <- stats[match(common, rownames(stats)),]
  stats <- cbind(stats, Full_ann)
  
  print("finding DMRs")
  
  #applying dmrff
  dmrs <- dmrff(estimate=stats$estimate,
                se=stats$se,
                p.value=stats$p.value,
                methylation=mval,
                chr=stats$chr,
                pos=stats$pos,
                maxgap=500,
                verbose=T)
  dmrs<-mutate(dmrs,fdr=p.adjust(p.value,method='fdr'))
  
  #We just keep regions with > 1 CpG site and fdr adjusted p < 0.05.
  dmrs <- dmrs[which(dmrs$fdr < 0.05 & dmrs$n >1),]
  
  if(nrow(dmrs)==0){stop("No DMRs associated with Factor at FDR<0.05, n CpG sites > 1")}
  
  sites <- dmrff.sites(dmrs, stats$chr, stats$pos)
  sites <- cbind(sites, stats[sites$site, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island","enhancer","OpenChromatin","TFBS", "estimate", "p.value")])
  sites <- cbind(sites, dmr=dmrs[sites$region,c("start","end","z","p.adjust")])
  
  #select probes that we can check expression in the array
  sites %>% filter(UCSC_RefGene_Name %in% rownames(RNA_file1))->select_probes
  #select genes with at least two CpG sites
  select_probes %>% group_by(UCSC_RefGene_Name) %>% summarise(n=length(UCSC_RefGene_Name)) %>% filter(n>1)->temp
  select_probes %>% filter(UCSC_RefGene_Name %in% temp$UCSC_RefGene_Name)->select_probes 
  
  print("Correlation between DMRs and Gene expression")
  
  results_corr=data.frame()
  for(i in unique(select_probes$UCSC_RefGene_Name)){
    #get mean values CpG sites within DMRs
    list_cpg=sites[which(sites$UCSC_RefGene_Name==i),5]
    temp_meth=data.frame(ID=colnames(betas),
                         avg_beta=colMeans(betas[rownames(betas) %in% list_cpg,],na.rm = T),
                         Feature=i)
    temp_meth=left_join(temp_meth,sampleSheet1 %>% dplyr::select(Factor,ID),by='ID')
    
    #get the expression for the same samples
    temp_exp=data.frame(ID=colnames(RNA_file1),
                        Expression=t(RNA_file1[rownames(RNA_file1)==i,]))
    names(temp_exp)[2]='Expression'
    
    #merge meth and exp
    temp_data=left_join(temp_meth,temp_exp)
    
    #run cor.test
    eval(parse(text=paste0("temp=cor.test(temp_data$avg_beta,temp_data$Expression)")))
    eval(parse(text=paste0("temp_r=cor.test(temp_data$avg_beta,temp_data$Factor)")))
    temp1=data.frame(Factor=MOFA_factors,feature=i,estimate_cor=temp$estimate,p_cor=temp$p.value,r2_Factor_Meth=temp_r$estimate*temp_r$estimate)
    results_corr=rbind(results_corr,temp1)
  }
    
    #processing
    results_corr %>% mutate(fdr_cor=p.adjust(p_cor,method='fdr'))->results_corr
    results_corr %>% filter(fdr_cor<0.05 & estimate_cor<0) ->temp#inverse correlations
    select_probes %>% mutate(Reg=ifelse(UCSC_RefGene_Name %in% temp$feature,'Yes','No'))->select_probes
    select_probes %>% mutate(CHR=extract_numeric(chr))->select_probes
    select_probes=left_join(select_probes,results_corr,by=c('UCSC_RefGene_Name'='feature'))
    
    #formatting cols
    select_probes %>% dplyr::select('Factor','region','chr','dmr.start','dmr.end',
                                    'estimate_cor','p_cor','fdr_cor',"UCSC_RefGene_Name","Name",
                                    "pos","UCSC_RefGene_Name","UCSC_RefGene_Group",
                                    "Relation_to_Island","enhancer","OpenChromatin",
                                    "TFBS","dmr.z","dmr.p.adjust")->select_probes
    #save results
    write_csv(select_probes, paste0(output_directory,MOFA_factors,"_DMR_Expression.csv"))
}
  
  
#output expected: CSV format file with MOFA Factor, DMR region, chromossome, start and end of DMRs, 
#                 correlation DMR and gene expression,p-value of correlation, fdr of correlation,
#                 Gene name, CpG site probe,hg19 position of each probe, followed by genomic annotations,
#                 dmr.z=z-score of DMR with MOFA factor,dmr.p.adjust=p-value of DMR with MOFA factor 






























