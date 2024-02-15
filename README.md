# Multi-omics analysis

R scripts used for the clear cell renal carcinoma multi-omics analyses.
R scripts using RStudio Server (2022.02.4) with R version 4.1.
The versions of each R package are reported at the beginning of code. 

## Main R scripts

### The R scripts for the core analyses are listed below

- [pre_processing_DNAMethyl.Rmd](MainAnalysis_code/pre_processing_DNAMethyl.Rmd): R code to process DNA methylation data. 
It was adapted from https://github.com/IARCbioinfo/methylkey.
1. It loads the idat files (probe intensities; unmethylated and methylated).
2. It uses minfi package to extract beta values.
3. It infers sex and age based on pre-selected CpG sites. 
4. It estimates sample-specific quality control (QC) for methylation data.
5. It normalises the beta values using different methods (e.g.,'Funmorm').
6. It filters out low-quality, highly missing, and blacklist-flagged CpG sites.
7. It corrects for potential batch effects (if 'Batch' variable present).
8. It outputs matrices containing beta and m-values, and pre-processed intensity files.


- [Format_Files_and_run_MOFA.R](MainAnalysis_code/Format_Files_and_run_MOFA.R): R code to format files and run MOFA.
1. It loads each Omic layer (DNAm, RNA, and WGS-somatic mutation profile).
2. It ranks genomic features by variance across samples within each Omic layer.
3. It selects the top 5,000 most-variable features of DNAm and RNA levels.
4. Plot the cumulative variance of the top features vs. all features within Omic layers.
5. Combine DNAm (DNA methylome), RNA (transcriptome), and WGS (WGS-based mutation signatures SBS, DBS, SV, CNV) into MOFA input as shown in the example from our study  ./Data/Example_MOFA_input.RData 
6. Run Multi-Omics Factor Analysis (MOFA)


- [LASSO_signatures.R](MainAnalysis_code/LASSO_signatures.R): R code to select the most informative features correlated with MOFA factors and used them to generate signatures to infer the MOFA factors in other datasets.
1. This function uses the MOFA object file (./Data/Example_MOFA_input.RData) as input, MOFA Factor (e.g., Factor1...Factor10), and Omic layer ('DNAm' for DNA methylome and 'RNA' for transcriptome).
2. Run Pearson's correlation test between genomic variables and MOFA Factors.
3. Select the MOFA Factor-correlated genomic features after multiple-testing correction (False Discovery rate <0.05).
4. It uses LASSO models (R tidymodels) to select the most informative features.
5. Outputs a data frame with LASSO coefficients per genomic variable that will be used to generate Factor signatures.


- [Inferring_MOFA_factors.R](MainAnalysis_code/Inferring_MOFA_factors.R): R code to infer MOFA Factor signatures in other datasets.
1. This function uses the output of LASSO_signatures.R for queried MOFA factor and path to DNA methylome (m-values) or RNA (log2-voom transformed if RNA-seq) matrices in TXT format, rows as feature, cols as samples.
2. Check the proportion of genomic features found overlapping with MOFA signatures (proceed if more than 75%).
3. Scale each genomic feature across samples (feature divided by standard deviation)
4. Multiply each genomic feature by respective LASSO coefficients.
5. Sum values of each genomic feature up by sample.
6. Calculate the variance in the queried Omic layer explained by MOFA factor signature.


- [DMR_Expression.R](MainAnalysis_code/DMR_Expression.R): R code to evaluate the relationship between differentially methylated regions and gene expression.
This function uses the outputs from pre_processing_DNAMethyl.Rmd (beta and m-values, RGSet) and a data frame with samples as rows and MOFA factor values with covariates as cols. This data frame should also include sequencing ID as 'Basename' and common id as 'ID'.
1. It loads DNA methylation and gene expression matrices (sample as cols, gene/CpG as rows) and replaces old sample id (Basename) by a common identifier to both DNA methylation and gene expression data. 
2. It creates an annotation file using minfi::getAnnotation for the CpG sites in RGSet.
3. It Filters a matrix of M-values (or beta values) by distance to SNP/variant. Also (optionally) removes cross-hybridising probes and sex-chromosome probes.
4. It runs linear regression models using limma (m-values ~ MOFA factor + covariates)# list of covariates could be adjusted if needed.
5. Select the DMRs with p.adjusted (FDR) values <0.05 with >1 CpG sites with tested region.
6. Average the betas within DMRs and evaluate the correlation with MOFA factor and gene expression levels.
7. It outputs an extensively annotated CSV file containing CpG sites per DMR, p-values and correlations between DMR, MOFA factor, and gene expression levels, and annotations of DMR.


- [Gobal_Methylation.R](MainAnalysis_code/Gobal_Methylation.R): This R code uses the output from pre_processing_DNAMethyl.Rmd to annotate and infer global methylation levels using Alu and Line-1 repetitive elements.
1. It loads RGSet and m-values.
2. Use REMP package to predict DNA methylation of locus-specific repetitive elements (RE) by learning surrounding genetic and epigenetic information.
3. It outputs a data frame in CSV format with mean and median methylation levels within Line-1 and Alu repetitive elements.


- [Pathway_Analyses_GSEA.R](MainAnalysis_code/Pathway_Analyses_GSEA.R): This function uses a data frame with the top correlated features of interest as input
1. Its reads annotation (GO and hallmarks) and correlation files. see examples in ./Data/
2. It runs the pathway analysis using fgsea::fgsea function.
3. It collapses redundant pathways and merge GO terms with hallmarks of cancer annotations.
4. It outputs a data frame in CSV format with the results of the pathway enrichment.













