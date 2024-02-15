###################################################
#calculating global methylation changes 
#BiocManager::install("REMP")

#this code used the RGSet file to annotate and m-values files (outputs from pre_processing_DNAMethyl.Rmd)

library(REMP)#v1.18.0
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)#v0.6.0
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)#v0.6.0

#load DNAm methylation file saved as rdata (RGSet, mval)
load(file = rda_path)

Full_ann=as.data.frame(minfi::getAnnotation(RGset, lociNames = NULL,
                                            orderByLocation = FALSE, dropNonMapping = FALSE))
#for Alu elements
Alu <- remprofile(mval, REtype = "Alu", annotation.source = "UCSC", genome = "hg19", Seq.GR = NULL, RE = NULL, impute = FALSE, imputebyrow = FALSE, verbose = FALSE)
details(Alu)
beta_Alu <- as.data.frame(rempB(Alu))
m_Alu <- as.data.frame(rempM(Alu))

#for LINE-1 elements
Line1 <- remprofile(mval, REtype = "L1", annotation.source = "UCSC", genome = "hg19", Seq.GR = NULL, RE = NULL, impute = FALSE, imputebyrow = FALSE, verbose = FALSE)
details(Line1)
beta_Line1 <- as.data.frame(rempB(Line1))
m_Line1 <- as.data.frame(rempM(Line1))

### Calculate mean Alu and Line1 methylation levels for each sample 
beta_Alu.Mean <- as.data.frame(colMeans(beta_Alu,na.rm=T))
beta_Alu.Mean$Basename <- rownames(beta_Alu.Mean)
colnames(beta_Alu.Mean)[1] <- "beta_Alu.Mean"

beta_Alu.Median <- as.data.frame(apply(beta_Alu,2,median))
beta_Alu.Median$Basename  <- rownames(beta_Alu.Median)
colnames(beta_Alu.Median)[1] <- "beta_Alu.Median"

m_Alu.Mean <- as.data.frame(colMeans(m_Alu,na.rm=T))
m_Alu.Mean$Basename  <- rownames(m_Alu.Mean)
colnames(m_Alu.Mean)[1] <- "m_Alu.Mean"

m_Alu.Median <- as.data.frame(apply(m_Alu,2,median))
m_Alu.Median$Basename  <- rownames(m_Alu.Median)
colnames(m_Alu.Median)[1] <- "m_Alu.Median"

beta_Line1.Mean <- as.data.frame(colMeans(beta_Line1,na.rm=T))
beta_Line1.Mean$Basename  <- rownames(beta_Line1.Mean)
colnames(beta_Line1.Mean)[1] <- "beta_Line1.Mean"

beta_Line1.Median <- as.data.frame(apply(beta_Line1,2,median))
beta_Line1.Median$Basename  <- rownames(beta_Line1.Median)
colnames(beta_Line1.Median)[1] <- "beta_Line1.Median"

m_Line1.Mean <- as.data.frame(colMeans(m_Line1,na.rm=T))
m_Line1.Mean$Basename  <- rownames(m_Line1.Mean)
colnames(m_Line1.Mean)[1] <- "m_Line1.Mean"

m_Line1.Median <- as.data.frame(apply(m_Line1,2,median))
m_Line1.Median$Basename  <- rownames(m_Line1.Median)
colnames(m_Line1.Median)[1] <- "m_Line1.Median"

df_list <- list(m_Alu.Mean,m_Alu.Median,m_Line1.Mean,m_Line1.Median,beta_Alu.Mean,beta_Alu.Median,beta_Line1.Mean,beta_Line1.Median)
global_meth <- Reduce(function(x,y) merge(x,y,by="Basename"), df_list)

global_meth$Basename=gsub('X','',global_meth$Basename)

#save results
write_csv(global_meth,'file_name.csv')
save(global_meth,Line1,Alu, file="file_name.RData")

#Expected output
#data frame with mean and median values of betas and m-values as cols, samples as rows.
