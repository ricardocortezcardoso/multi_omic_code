#DMR using DNA methylome and how it relates to gene expression
#BiocManager::install("DMRcate")
#devtools::install_github("perishky/dmrff")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

##DNAm_file: path to rda file with processed (betas and mval) and and RGset methylation files
##RNA_file: path to rd file with processed gene expression data (genes as rows, sample as cols)
##samplesheet_file: info file wit patients data (features as cols, samples as rows), including sequencing (Basename) and epi names (ID)


#load packages
library(tidyverse)
library(data.table)
library(dmrff)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(limma)


calculate_DMR <- function(DNAm_file, RNA_file, samplesheet_file, MOFA_factor, output_directory) {
  
  #load normalized beta files
  load(DNAm_file)
  betas=as.data.frame(betas)
  data.table::setnames(betas, old=samplesheet_file$Basename, new=samplesheet_file$ID,skip_absent=TRUE)
  betas=as.matrix(betas)
  
  mval=as.data.frame(mval)
  data.table::setnames(mval, old=samplesheet_file$Basename, new=samplesheet_file$ID,skip_absent=TRUE)
  mval=as.matrix(mval)
  
  #get annotation file
  Full_ann=as.data.frame(minfi::getAnnotation(RGset, lociNames = NULL,
                                              orderByLocation = FALSE, dropNonMapping = FALSE))
  Full_ann %>%  
    mutate(UCSC_RefGene_Name=gsub(";.*",'',UCSC_RefGene_Name),
           UCSC_RefGene_Group=gsub(";.*",'',UCSC_RefGene_Group),
           enhancer=ifelse(nchar(Phantom4_Enhancers)>1 | nchar(Phantom5_Enhancers)>1 | nchar(X450k_Enhancer)>1 ,'Yes','No'),
           OpenChromatin=ifelse(!OpenChromatin_Evidence_Count=='','Open_Chromatin',''),
           TFBS=ifelse(!TFBS_Evidence_Count=='','TFBS',''))->Full_ann
  
  #load RNA file
  load(RNA_file)
  RNA_file=RNA_file[,colnames(RNA_file) %in% colnames(mval)]
  
  #final QC for methylation data
  noSNPs <- DMRcate::rmSNPandCH(betas, dist=2, mafcut=0.05)
  nrow(noSNPs)
  mval <- mval[rownames(noSNPs),] 
  colnames(noSNPs) <- colnames(mval)
  
  #running limma
  design <- model.matrix(~ MOFA_factor+age_diag+sex+country,data = samplesheet_file)
  fit <- limma::lmFit(mval, design)
  fit <- limma::eBayes(fit)
  stats<-data.frame(estimate=fit$coefficients,
                    se=sqrt(fit$s2.post) * fit$stdev.unscaled,
                    p.value=fit$p.value)
  
  #creating annotation file
  common <- intersect(rownames(mval), rownames(Full_ann))
  Full_ann <- Full_ann[match(common, rownames(Full_ann)),]
  stats <- stats[match(common, rownames(stats)),]
  stats <- cbind(stats, Full_ann)
  
  ###applying dmrff
  dmrs <- dmrff(estimate=stats$estimate.MOFA_factor,
                se=stats$se.MOFA_factor,
                p.value=stats$p.value.MOFA_factor,
                methylation=mval,
                chr=stats$chr,
                pos=stats$pos,
                maxgap=500,
                verbose=T)
  
  #We just keep regions with > 1 CpG site and Bonferroni adjusted p < 0.05.
  dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]
  sites <- dmrff.sites(dmrs, stats$chr, stats$pos)
  sites <- cbind(sites, stats[sites$site, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island","enhancer","OpenChromatin","TFBS", paste0("estimate.",MOFA_factor), paste0("p.value.",MOFA_factor))])
  sites <- cbind(sites, dmr=dmrs[sites$region,c("start","end","z","p.adjust")])
  
  #select probes that we can check expression in the array
  sites %>% filter(UCSC_RefGene_Name %in% rownames(RNA_file))->select_probes
  #select genes with at least two CpG sites
  select_probes %>% group_by(UCSC_RefGene_Name) %>% summarise(n=length(UCSC_RefGene_Name)) %>% filter(n>1)->temp
  select_probes %>% filter(UCSC_RefGene_Name %in% temp$UCSC_RefGene_Name)->select_probes 
  
  results_corr=data.frame()
  for(i in unique(select_probes$UCSC_RefGene_Name)){
    #get mean values CpG sites within DMRs
    list_cpg=sites[which(sites$UCSC_RefGene_Name==i),5]
    temp_meth=data.frame(ID=colnames(betas),
                         avg_beta=colMeans(betas[rownames(betas) %in% list_cpg,],na.rm = T),
                         Feature=i)
    temp_meth=left_join(temp_meth,samplesheet_file %>% dplyr::select(MOFA_factor,ID),by='ID')
    
    #get the expression for the same samples
    temp_exp=data.frame(ID=colnames(RNA_file),
                        Expression=t(RNA_file[rownames(RNA_file)==i,]))
    names(temp_exp)[2]='Expression'
    
    #merge meth and exp
    temp_data=left_join(temp_meth,temp_exp)
    
    #run cor.test
    eval(parse(text=paste0("temp=cor.test(temp_data$avg_beta,temp_data$Expression)")))
    eval(parse(text=paste0("temp_r=cor.test(temp_data$avg_beta,temp_data$",MOFA_factor,")")))
    temp1=data.frame(Factor=MOFA_factor,feature=i,estimate=temp$estimate,p=temp$p.value,r2_Factor_Meth=temp_r$estimate*temp_r$estimate)
    results_corr=rbind(results_corr,temp1)
  }
    
    #processing
    results_corr %>% mutate(fdr=p.adjust(p,method='fdr'))->results_corr
    results_corr %>% filter(fdr<0.05 & estimate<0) ->temp
    select_probes %>% mutate(Reg=ifelse(UCSC_RefGene_Name %in% temp$feature,'Yes','No'))->select_probes
    select_probes %>% mutate(CHR=extract_numeric(chr))->select_probes
    select_probes=left_join(select_probes,results_corr,by=c('UCSC_RefGene_Name'='feature'))
    
    #save results
    write_csv(select_probes, paste0(output_directory,MOFA_factor,"_DMR_Expression.csv"))
}
  
  
  
  
  
  
 
 















###load normalized beta files
load("/data/geniluc/work/cortezr/RCC_Mutographs/RCC_121_Funnorm_2023-09-04.rda")
betas=as.data.frame(betas)
data.table::setnames(betas, old=sampleSheet$Basename, new=sampleSheet$Sample_ID,skip_absent=TRUE)
betas=as.matrix(betas)

mval=as.data.frame(mval)
data.table::setnames(mval, old=sampleSheet$Basename, new=sampleSheet$Sample_ID,skip_absent=TRUE)
mval=as.matrix(mval)

#get annotation file
Full_ann=as.data.frame(minfi::getAnnotation(RGset, lociNames = NULL,
                                            orderByLocation = FALSE, dropNonMapping = FALSE))
Full_ann %>%  
  mutate(UCSC_RefGene_Name=gsub(";.*",'',UCSC_RefGene_Name),
         UCSC_RefGene_Group=gsub(";.*",'',UCSC_RefGene_Group),
         enhancer=ifelse(nchar(Phantom4_Enhancers)>1 | nchar(Phantom5_Enhancers)>1 | nchar(X450k_Enhancer)>1 ,'Yes','No'),
         OpenChromatin=ifelse(!OpenChromatin_Evidence_Count=='','Open_Chromatin',''),
         TFBS=ifelse(!TFBS_Evidence_Count=='','TFBS',''))->Full_ann

#load IARC transcriptome data
load('IARC_Expression151RCC_QC.rda')
IARC_expression_profile=IARC_expression_profile[,colnames(IARC_expression_profile) %in% colnames(mval)]

#final QC for methylation data
noSNPs <- DMRcate::rmSNPandCH(betas, dist=2, mafcut=0.05)#5836 probes removed  because of 2 bases near SNPs
nrow(noSNPs)
mval <- mval[rownames(noSNPs),] #filter Mset by the set of probes to be analysed. Probes already normalized and QC 
colnames(noSNPs) <- colnames(mval)

#get the variable I want to contrast: Factor1,2,6
LF_Discovery<-as.data.frame(read_csv("DiscoverySet_04102023.csv"))
samplesheet_RCC=data.frame(other_id=colnames(mval))
samplesheet_RCC=left_join(samplesheet_RCC,LF_Discovery %>% dplyr::select(other_id,Factor1,Factor2,Factor6,age_diag,sex,country))

all(samplesheet_RCC$other_id==colnames(mval))#TRUE

#running limma
design <- model.matrix(~Factor6+age_diag+sex+country,data = samplesheet_RCC)
fit <- limma::lmFit(mval, design)
fit <- limma::eBayes(fit)
stats<-data.frame(estimate=fit$coefficients,
                  se=sqrt(fit$s2.post) * fit$stdev.unscaled,
                  p.value=fit$p.value)


#creating annotation file
common <- intersect(rownames(mval), rownames(Full_ann))
Full_ann <- Full_ann[match(common, rownames(Full_ann)),]
stats <- stats[match(common, rownames(stats)),]
stats <- cbind(stats, Full_ann)

###applying dmrff
dmrs <- dmrff(estimate=stats$estimate.Factor6,
              se=stats$se.Factor6,
              p.value=stats$p.value.Factor6,
              methylation=mval,
              chr=stats$chr,
              pos=stats$pos,
              maxgap=500,
              verbose=T)
#We just keep regions with > 1 CpG site and Bonferroni adjusted p < 0.05.
dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]
sites <- dmrff.sites(dmrs, stats$chr, stats$pos)
sites <- cbind(sites, stats[sites$site, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island","enhancer","OpenChromatin","TFBS","estimate.Factor6", "p.value.Factor6")])
sites <- cbind(sites, dmr=dmrs[sites$region,c("start","end","z","p.adjust")])

#select probes that we can check expression in the array
sites %>% filter(UCSC_RefGene_Name %in% rownames(IARC_expression_profile))->select_probes
#select genes with at least two CpG sites
select_probes %>% group_by(UCSC_RefGene_Name) %>% summarise(n=length(UCSC_RefGene_Name)) %>% filter(n>1)->temp
select_probes %>% filter(UCSC_RefGene_Name %in% temp$UCSC_RefGene_Name)->select_probes

results_corr=data.frame()
for(i in unique(select_probes$UCSC_RefGene_Name)){
  #get mean values CpG sites within DMRs
  list_cpg=sites[which(sites$UCSC_RefGene_Name==i),5]
  temp_meth=data.frame(IARC_ID=colnames(betas),
                       avg_beta=colMeans(betas[rownames(betas) %in% list_cpg,],na.rm = T),
                       Feature=i)
  temp_meth=left_join(temp_meth,LF_Discovery %>% dplyr::select(Factor6,other_id),by=c('IARC_ID'='other_id'))
  
  #get the expression for the same samples
  temp_exp=data.frame(IARC_ID=colnames(IARC_expression_profile),
                      Expression=t(IARC_expression_profile[rownames(IARC_expression_profile)==i,]))
  names(temp_exp)[2]='Expression'
  
  #merge meth and exp
  temp_data=left_join(temp_meth,temp_exp)
  
  #run cor.test
  eval(parse(text=paste0("temp=cor.test(temp_data$avg_beta,temp_data$Expression)")))
  eval(parse(text=paste0("temp_r=cor.test(temp_data$avg_beta,temp_data$Factor6)")))
  temp1=data.frame(Factor='Factor6',feature=i,estimate=temp$estimate,p=temp$p.value,r2_Factor_Meth=temp_r$estimate*temp_r$estimate)
  results_corr=rbind(results_corr,temp1)
}

results_corr %>% mutate(fdr=p.adjust(p,method='fdr'))->results_corr
results_corr %>% filter(fdr<0.05 & estimate<0) ->temp
select_probes %>% mutate(Reg=ifelse(UCSC_RefGene_Name %in% temp$feature,'Yes','No'))->select_probes
select_probes %>% mutate(CHR=extract_numeric(chr))->select_probes
select_probes=left_join(select_probes,results_corr,by=c('UCSC_RefGene_Name'='feature'))
write_csv(select_probes,'./Results/DMRsFactor6_Expression.csv')

