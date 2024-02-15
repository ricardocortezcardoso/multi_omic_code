# Processing File for MOFA

# Load libraries
library(tidyverse)#v2.0.0
library(data.table)#v1.14.8
library(bigreadr)#v0.2.5
library(MOFA2)#v1.9.2

## DNA Methylation Files

# Import files with the normalized m-values for DNA methylome (rda_path)
load(file = rda_path)

# Remove NAs introduced to replace low-quality sites
probes <- names(which(rowSums(apply(betas, 2, is.na)) > 0))
betas <- betas[!rownames(betas) %in% probes,]

# Filter m-table (m_filt) for minimum beta variance (>0.1)
beta.diff <- apply(betas, 1, max, na.rm = TRUE) - apply(betas, 1, min, na.rm = TRUE)
names(beta.diff) <- rownames(betas)
beta.diff <- beta.diff[beta.diff > 0.1]

beta_filt <- betas[which(rownames(betas) %in% names(beta.diff)), ]
dim(beta_filt)
m_filt <- mval[which(rownames(mval1) %in% names(beta.diff)), ]

# Check the correct probes have been removed
beta_filt.diff <- apply(beta_filt, 1, max, na.rm = TRUE) - apply(beta_filt, 1, min, na.rm = TRUE)
plot(density(beta_filt.diff))  # zero density at 0.1

# Reduce m-table to top 5,000 most variable probes
vv <- apply(m_filt, 1, var)
cv <- cumsum(sort(vv, decreasing = TRUE)) / sum(vv)
all(order(cv, decreasing = FALSE) == 1:length(cv))

m_filt_top <- m_filt[order(match(rownames(m_filt), names(cv))), ]
m_filt_top <- m_filt_top[1:5000, ]
beta_filt_top <- beta_filt[rownames(beta_filt) %in% rownames(m_filt_top), ]

# Plot M-values of top 5,000 CpG sites
png(paste0(output_directory,'/Distribution_topMeth.png'), units = "in", width = 4, height = 4, res = 600)
minfi::densityPlot(m_filt_top, xlab = 'm-values', pal = '#BC3C29FF', main = 'M-values of top 5,000 CpG sites')
dev.off()

# Plotting the Variance of all CpGs (EPIC) versus top 5,000 selected for MOFA
temp <- data.frame(Variance = cv)
temp$feature <- rownames(temp)
temp %>% mutate(Group = ifelse(feature %in% rownames(beta_filt_top), 'Yes', 'No')) -> temp
temp$feature <- fct_reorder(temp$feature, temp$Variance)
temp$Variance[5000]

to_plot <- ggplot(data = temp, aes(x = feature, y = Variance * 100, color = Group)) +
  scale_color_manual(values = c('#808180FF', '#BC3C29FF')) +
  geom_jitter(size = 2, na.rm = TRUE) +
  geom_hline(yintercept = temp$Variance[5000] * 100, linetype = 'dashed', color = 'black') +
  xlab('CpG sites (N=707,610)') +
  ylab('Cumulative Variance of CpG sites (%)') +
  scale_x_discrete(drop = FALSE, breaks = seq(0, 707610, 1000), na.translate = TRUE) +
  scale_y_continuous(breaks = seq(0, 100, 5), limits = c(0, 100)) +
  ggtitle('Cumulative Variance of CpG sites in 120 ccRCC tumors') +
  theme_classic(base_size = 30) +
  theme(legend.position = 'top', axis.text.x = element_blank(), plot.title = element_text(size = 28, face = 'bold')) +
  labs(color = 'Top 5,000 CpGs')

ggsave(plot = to_plot, filename = paste0(output_directory,"CV_CpGs_120ccRCC.png"), device = "png", width = 15, height = 12, dpi = 600)



## Transcriptome file (Probe intensity of array as gene expression)

# Open expression profile files (transcriptome): probe intensity (RNA_data: probes as rows and samples as cols), patient info (RNA_info), and probe info (RNA_probe)
RNA_data <- bigreadr::fread2(path_to_RNA_data)
RNA_probe <- bigreadr::fread2(path_to_RNA_probe)#probe.qc (detection rate),Probe_Chr(chr),Gene_Symbol (Gene Symbol)
RNA_info <- bigreadr::fread2(path_to_RNA_info)

# Keep good quality probes: call detection rate > 5% in both paired Cases and Controls, choosing the ones with the highest qc per gene
good_probes<-RNA_probe %>%
  filter(probe.qc > 0.05 & Probe_Chr %in% c(paste0('chr', seq(1,22)))) %>%
  group_by(Gene_Symbol) %>%
  arrange(desc(probe.qc)) %>%
  dplyr::slice(1)


RNA_data <- RNA_data[rownames(RNA_data) %in% good_probes$Probe_ID, ]
RNA_data <- na.omit(RNA_data)
RNA_data <- RNA_data[good_probes$Probe_ID, ]  # Probes same order as good_probes file
all(good_probes$Probe_ID == rownames(RNA_data))  # TRUE
rownames(IARC_expression_profile) <- good_probes$Gene_Symbol

# Filter expression for minimum FPKM variance (>1 FPKM)
expression_diff <- apply(RNA_data, 1, max, na.rm = TRUE) - apply(RNA_data, 1, min, na.rm = TRUE)
names(expression_diff) <- rownames(RNA_data)
expression_diff <- expression_diff[which(expression_diff >= 1)]  # 11,545 genes
expression_filt <- RNA_data[rownames(RNA_data) %in% names(expression_diff), ]
dim(expression_filt)

# Reduce expression_filt to top 5,000 most variable genes across samples
vv <- apply(expression_filt, 1, var)
cv <- cumsum(sort(vv, decreasing = TRUE)) / sum(vv)
all(order(cv, decreasing = FALSE) == 1:length(cv))  # TRUE
expression_filt <- expression_filt[order(match(rownames(expression_filt), names(cv))), ]
expression_filt_top <- expression_filt[1:5000, ]

expression_filt_top <- as.matrix(expression_filt_top)

# Plotting the Variance of all gene expression versus top 5,000 selected for MOFA
temp <- data.frame(Variance = cv)
temp$feature <- rownames(temp)
temp %>% mutate(Group = ifelse(feature %in% rownames(expression_filt_top), 'Yes', 'No')) -> temp
temp$feature <- fct_reorder(temp$feature, temp$Variance)

to_plot <- ggplot(data = temp, aes(x = feature, y = Variance * 100, color = Group)) +
  scale_color_manual(values = c('#808180FF', '#BC3C29FF')) +
  geom_jitter(size = 2, na.rm = TRUE) +
  geom_hline(yintercept = temp$Variance[5000] * 100, linetype = 'dashed', color = 'black') +
  xlab('RNA levels (N=11,545)') +
  ylab('Cumulative Variance of RNA levels (%)') +
  scale_y_continuous(breaks = seq(0, 100, 5), limits = c(0, 100)) +
  ggtitle('Cumulative Variance of RNA levels in 151 ccRCC tumors') +
  theme_classic(base_size = 28) +
  theme(legend.position = 'top', axis.text.x = element_blank(), plot.title = element_text(size = 28, face = 'bold')) +
  labs(color = 'Top 5,000 Genes')

ggsave(plot = to_plot, filename = paste0(output_directory,"CV_RNA_151ccRCC.png"), device = "png", width = 15, height = 12, dpi = 600)

png(paste0(output_directory,"Distribution_topRNA.png"), units = "in", width = 4, height = 4, res = 600)
minfi::densityPlot(expression_filt_top, xlab = 'RNA levels', pal = '#BC3C29FF', main = 'Top 5,000 RNA levels')
dev.off()



##Somatic profile

#open somatic file (features as cols, samples as rows) with driver mutations (0 absence, 1 present), DNA mutational signatures (continuous) 
somatic_mut <- as.data.frame(read_csv(path_to_somatic_file))

#counting the non-zeros, and exclude features < 5 events across samples
nonzero <- function(x) sum(x != 0, na.rm = TRUE)
to_exc=as.data.frame(t(somatic_mut %>% summarise_all(nonzero)))
to_exc %>% filter(V1<5)->to_exc 
somatic_mut %>% dplyr::select(!rownames(to_exc))->somatic_mut

WGS=as.matrix(t(somatic_mut))



##Match sample names across Omic layers
Combined_samples <- sort(unique(c(colnames(expression_filt_top),colnames(m_filt_top),colnames(WGS))))

# Format transcriptome, DNA methylation data and somatic mutation profile as MOFA inputs
D_expr_MOFA <- matrix(NA,nrow(expression_filt_top),length(which(!Combined_samples%in%colnames(expression_filt_top))),dimnames = list(rownames(expression_filt_top),Combined_samples[which(!Combined_samples%in%colnames(expression_filt_top))]))
D_expr_MOFA <- as.data.frame(D_expr_MOFA) 
D_expr_MOFA <- merge(expression_filt_top,D_expr_MOFA,by.x="row.names", by.y="row.names")
D_expr_MOFA <- data.frame(D_expr_MOFA[,-1],row.names=D_expr_MOFA[,1])
D_expr_MOFA <- D_expr_MOFA[,order(match(colnames(D_expr_MOFA),Combined_samples))]
dim(D_expr_MOFA) 

D_met_MOFA <- matrix(NA, nrow(m_filt_top), length(which(!Combined_samples %in% colnames(m_filt_top))), dimnames = list(rownames(m_filt_top), Combined_samples[which(!Combined_samples %in% colnames(m_filt_top))]))
D_met_MOFA <- as.data.frame(D_met_MOFA) 
D_met_MOFA <- merge(m_filt_top, D_met_MOFA, by.x="row.names", by.y="row.names")
D_met_MOFA <- data.frame(D_met_MOFA[,-1], row.names=D_met_MOFA[,1])
D_met_MOFA <- D_met_MOFA[, order(match(colnames(D_met_MOFA), Combined_samples))]
dim(D_met_MOFA) 

D_WGS_MOFA <- WGS[, order(match(colnames(WGS), Combined_samples))]
dim(D_WGS_MOFA) 

#quick check
all(colnames(D_expr_MOFA)==Combined_samples)
all(colnames(D_met_MOFA)==Combined_samples)
all(colnames(D_WGS_MOFA)==Combined_samples)

#convert to matrices format
D_expr_MOFA=as.matrix(D_expr_MOFA)
D_met_MOFA=as.matrix(D_met_MOFA)
D_WGS_MOFA=as.matrix(D_WGS_MOFA)


##Run MOFA

#settings
MOFAobject_ABC <- MOFA2::create_mofa(list("RNA" = D_expr_MOFA, "DNAm" = D_met_MOFA, "WGS" = D_WGS_MOFA))#example provided
data_optsB <- MOFA2::get_default_data_options(MOFAobject_ABC)
model_optsB <- MOFA2::get_default_model_options(MOFAobject_ABC)
model_optsB$num_factors <- 10
train_opts <- MOFA2::get_default_training_options(MOFAobject_ABC)
train_opts$convergence_mode <- "slow"
stochastic_opts <- MOFA2::get_default_stochastic_options(MOFAobject_ABC)

# Prepare
MOFAobject_ABC <- MOFA2::prepare_mofa(object=MOFAobject_ABC, data_options=data_optsB, model_options=model_optsB, training_options=train_opts)
save(MOFAobject_ABC, file = paste0(output_directory,"Example_MOFA_input.RData"))

#run
MOFAobject.trained <- MOFA2::run_mofa(MOFAobject_ABC, save_data = T, outfile = paste0(output_directory,"MOFAobject_trained.hdf5"),use_basilisk = T)

