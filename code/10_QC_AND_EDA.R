########################################################################################################################
#
#   10_QC_AND_EDA.R 
#
#   Perform initial quality control and exploratory data analysis on 4 y/o non-human primate RNA-seq samples
#   from DLPFC, ACC, HC, and V1
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(corrplot)
library(ggplot2)
library(WGCNA)
library(edgeR)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

pdf('./docs/10_QC_AND_EDA_pt1.pdf')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Retrieve datExpr, datMeta, and clean dataframes
#   (2) Visualize and filter raw read counts
#   (3) Normalize data and visualize trait correlations
#   (4) Outlier removal by connectivity z-score
#   (5) Export dataframes to RData files
#
########################################################################################################################

### (1) Retrieve datExpr, datMeta, and clean dataframes ################################################################

# raw count data from txImport provided by M Gandal with pseudo-alignment to macaque genome using Kallisto
load('./raw_data/txi.kallisto.MMUL8v87.rda')
mmul_annotEns87 <- annotEns87
mmul_txi.gene <- txi.gene
mmul_txi.tx <- txi.tx

# raw count data from txImport provided by M Gandal with pseudo-alignment to human genome using Kallisto
load('./raw_data/txi.kallisto.Hg38v84.rda')
hg_annotEns84 <- annotEns84
hg_txi.gene <- txi.gene
hg_txi.tx <- txi.tx

# import sample meta data and clean dataframe
datMeta <- read.csv(file = './raw_data/datMeta_20160601.csv', header = TRUE, sep = ',')
rownames(datMeta) <- datMeta$Sample
datMeta <- datMeta[,c(5:11,16:21)]
colnames(datMeta) <- c(colnames(datMeta)[1:8], 'MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4', 'MG_SeqDepth')
# for the seqPCs, MG stands for M Gandal, since he calculated the seqPCs differently than NP

# import sequencing statistics from picard tools to calculate seqPCs
GCsummary <- read.csv(file = './raw_data/qc_picard/GCsummary.csv', header = TRUE, sep = ',')
PicardToolsQC <- read.csv(file = './raw_data/qc_picard/PicardToolsQC.csv', header = TRUE, sep = ',')
RNAseqAlign <- read.csv(file = './raw_data/qc_picard/RNAseqAlign.csv', header = TRUE, sep = ',')
RNAseqDuplication <- read.csv(file = './raw_data/qc_picard/RNAseqDuplication.csv', header = TRUE, sep = ',')
RNAseqQC <- read.csv(file = './raw_data/qc_picard/RNAseqQC.csv', header = TRUE, sep = ',')

# combine all sequencing statistics into a single dataframe
combined_seqStats <- cbind(GCsummary, PicardToolsQC, RNAseqAlign, RNAseqDuplication, RNAseqQC)
combined_seqStats <- combined_seqStats[, !is.na(colSums(combined_seqStats != 0))] # Remove columns with all NA
combined_seqStats <- combined_seqStats[, colSums(combined_seqStats != 0) > 0] # Remove columns where all values are zero
combined_seqStats <- combined_seqStats[!duplicated(lapply(combined_seqStats, summary))] # Remove duplicated columns

# remove seqStats that are highly correlated with eachother before computing seqPCs
# NOTE: the reason we want to remove highly correlated seqPCs in beacuse the increase the 
# contribution of their common underlying factor to the PCA (See link below)
# https://stats.stackexchange.com/questions/50537/should-one-remove-highly-correlated-variables-before-doing-pca
seqStat_correlation <- cor(combined_seqStats[2:39])
corrplot(seqStat_correlation, tl.cex = 0.5) # before filtering highly correlated sequencing statistics
rownames(combined_seqStats) <- gsub('.{20}$', '', combined_seqStats$X)
combined_seqStats <- combined_seqStats[, c(5:7,17,23,25,28,31:32,34:36,38)]
seqStat_correlation <- cor(combined_seqStats[1:13])
corrplot(seqStat_correlation) # after removing highly correlated sequencing statistics

# calculate seqPCs after removing highly correlated sequencing metadata
seqPCs <- prcomp(t(scale(combined_seqStats, scale = TRUE)), center = FALSE)
summary(seqPCs)
eigs <- seqPCs$sdev ^ 2 # compute the variance explained by each PC
seqPCs_var_explained <- round((eigs / sum(eigs)) * 100, 1) # compute the variance explained by each PC
seqPCs <- data.frame(seqPCs$rotation)

# plot the first two sequencing principle components
ggplot(seqPCs$rotation, aes(x = seqPCs$PC1, y = seqPCs$PC2, col = as.factor(datMeta$SeqBatch))) +
  ggtitle('\nSequencing Principle Components\n') + 
  geom_point(size = 5) +
  xlab(paste0('\nSeqPC1 (', seqPCs_var_explained[1], '%)\n' )) +
  ylab(paste0('\nSeqPC2 (', seqPCs_var_explained[2], '%)\n' )) +
  theme(plot.title  = element_text(hjust=.5, size = 24), axis.title=element_text(size=14)) + 
  labs(color = 'SeqBatch')

### Double checks that the order of the samples is the same as the metadata
### THIS STEP IS VERY IMPORTANT!!!
print(identical(colnames(hg_txi.gene$counts), colnames(mmul_txi.gene$counts)))
print(identical(colnames(mmul_txi.gene$counts), rownames(datMeta)))
print(identical(rownames(datMeta), rownames(seqPCs)))

# add new seqPCs to metadata label NP, these are the seqPCs calculated by N Page and create factors
colnames(seqPCs) <- paste0('NP_seq', colnames(seqPCs))
datMeta <- cbind(datMeta, seqPCs[,c(1:5)])
datMeta$Subject <- as.factor(datMeta$Subject)
datMeta$Region <- factor(datMeta$Region, levels = c('DLPFC', 'ACC', 'HC', 'V1'))
datMeta$Group <- factor(datMeta$Group, levels = c('con', 'poly1', 'poly2'))
datMeta$SeqBatch <- factor(datMeta$SeqBatch, levels = c('1', '2'))
datMeta$MIA <- factor(datMeta$MIA, levels = c('CON', 'MIA'))


### (2) Visualize and filter raw read counts #########################################################################

# creates lists of datExpr dataframes so that they can be looped through with the mapply function below
datExpr_list <- list(mmul_txi.gene$counts, mmul_txi.tx$counts, hg_txi.gene$counts, hg_txi.tx$counts)
datExpr_name_list <- list('mmul_gene raw_count density pre-normalizaton', 'mmul_tx raw_count density pre-normalizaton', 'hg_gene raw_count density pre-normalizaton', 'hg_tx raw_count density pre-normalizaton')

# loop through each datExpr dataframe and plot the raw_count density across each sample colored by sequencing batch
mapply(function(datExpr, datExpr_name) {
  log2_datExpr <- log2(datExpr + 1)
  i=1; plot(density((log2_datExpr),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
            main = datExpr_name, xlab="log2(raw_counts + 1)", xlim=c(-1,20), ylim=c(0,0.3))
  for(i in 2:dim(log2_datExpr)[2]){
    lines(density((log2_datExpr[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
  }
  legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = rev(as.numeric(as.factor(levels(datMeta$SeqBatch)))))
}, datExpr_list, datExpr_name_list)


# Identify and filter low expressed genes --> These tend to be very noisy and decrease power for overall DGE by increase multiple comparisons testing burden
# There are many ways to identify 'low expressed' genes; e.g., TPM/CPM > 0.1 in 50% of samples or counts > 10 in 25% samples, etc (a bit arbitrary)
genes_to_keep <- filterByExpr(DGEList(mmul_txi.gene$counts)) #need to mention this in the methods section
table(genes_to_keep)
mmul_txi.gene <- mmul_txi.gene$counts[genes_to_keep,]

genes_to_keep <- filterByExpr(DGEList(mmul_txi.tx$counts)) #need to mention this in the methods section
table(genes_to_keep)
mmul_txi.tx <- mmul_txi.tx$counts[genes_to_keep,]

genes_to_keep <- filterByExpr(DGEList(hg_txi.gene$counts)) #need to mention this in the methods section
table(genes_to_keep)
hg_txi.gene <- hg_txi.gene$counts[genes_to_keep,]

genes_to_keep <- filterByExpr(DGEList(hg_txi.tx$counts)) #need to mention this in the methods section
table(genes_to_keep)
hg_txi.tx <- hg_txi.tx$counts[genes_to_keep,]


# creates lists of filtered datExpr dataframes so that they can be looped through with the mapply function below
datExpr_list <- list(mmul_txi.gene, mmul_txi.tx, hg_txi.gene, hg_txi.tx)
datExpr_name_list <- list('mmul_gene filtered raw_count density pre-normalizaton', 'mmul_tx filtered raw_count density pre-normalizaton', 'hg_gene filtered raw_count density pre-normalizaton', 'hg_tx filtered raw_count density pre-normalizaton')

# loop through each datExpr dataframe and plot the raw_count density across each sample colored by sequencing batch
mapply(function(datExpr, datExpr_name) {
  log2_datExpr <- log2(datExpr + 1)
  i=1; plot(density((log2_datExpr),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
            main = datExpr_name, xlab="log2(raw_counts + 1)", xlim=c(-1,20), ylim=c(0,0.3))
  for(i in 2:dim(log2_datExpr)[2]){
    lines(density((log2_datExpr[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
  }
  legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = rev(as.numeric(as.factor(levels(datMeta$SeqBatch)))))
}, datExpr_list, datExpr_name_list)


### (3) Normalize data and visualize trait correlations ##############################################################

# perform TMM normalization on all datasets
mmul_txi.gene.norm <- cpm(calcNormFactors(DGEList(mmul_txi.gene), method = 'TMM'), log = TRUE)
mmul_txi.tx.norm <- cpm(calcNormFactors(DGEList(mmul_txi.tx), method = 'TMM'), log = TRUE)
hg_txi.gene.norm <- cpm(calcNormFactors(DGEList(hg_txi.gene), method = 'TMM'), log = TRUE)
hg_txi.tx.norm <- cpm(calcNormFactors(DGEList(hg_txi.tx), method = 'TMM'), log = TRUE)

# creates lists of filtered datExpr dataframes so that they can be looped through with the mapply function below
datExpr.norm_list <- list(mmul_txi.gene.norm, mmul_txi.tx.norm, hg_txi.gene.norm, hg_txi.tx.norm)
datExpr.norm_name_list <- list('mmul_gene_post-normalizaton', 'mmul_tx_post-normalizaton', 'hg_gene_post-normalizaton', 'hg_tx_post-normalizaton')

# loop through each datExpr dataframe and plot the raw_count density across each sample colored by sequencing batch
mapply(function(datExpr, datExpr_name) {
  i=1; plot(density((datExpr),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
            main = datExpr_name, xlab="cpm((TMM normalized counts), log = TRUE)", xlim=c(-5,20), ylim=c(0,0.20))
  for(i in 2:dim(datExpr)[2]){
    lines(density((datExpr[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
  }
  legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = rev(as.numeric(as.factor(levels(datMeta$SeqBatch)))))
}, datExpr.norm_list, datExpr.norm_name_list)


# Correlation panel function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if ((class(x) == "numeric" | class(x) == "integer") & (class(y) == "numeric" | class(y) == "integer")) {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.2/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


# loop through each expression dataframe and generate trait correlation plots
mapply(function(datExpr.norm, datExpr.norm_name) {
  
  pairsDat <- datMeta[,c(2:4,6:18)]
  figureDir <- './docs/trait_correlations/'
  
  ePCs <- prcomp(t(scale(t(datExpr.norm), scale = TRUE)), center = FALSE)
  ePCs <- data.frame(ePCs$rotation)

  pdf(paste0(figureDir, datExpr.norm_name, '_trait_pearson_correlations.pdf'), 24, 20)
  pairs(cbind(pairsDat, ePCs[,c(1:5)]), pch = 19, upper.panel = panel.cor, main = paste("Covariates and MaxQuant Comparison -- |Pearson's rho| correlation values","\n", datExpr.norm_name, sep = ""),cex.labels=1.2)
  dev.off()
  
}, datExpr.norm_list, datExpr.norm_name_list)


pdf('./docs/10_QC_AND_EDA_pt2.pdf')

# loop through each datExpr dataframe and plot the expression PCs by sequencing batch
mapply(function(datExpr.norm, datExpr.norm_name) {
  
  mds = cmdscale(dist(t(datExpr.norm)),k = 10)
  colnames(mds) = paste0("PC",1:ncol(mds))
  
  pairs(mds, col=factor(datMeta$SeqBatch), main=paste0('Expression PCs of ', datExpr.norm_name), pch=19)
  par(xpd = TRUE,oma=c(1,1,1,1)); legend('bottomright', levels(factor(datMeta$SeqBatch)),fill=1:2,cex=.5)
  
}, datExpr.norm_list, datExpr.norm_name_list)


### (4) Outlier removal by connectivity z-score ######################################################################

# creates lists of filtered datExpr dataframes so that they can be looped through with the mapply function below
datExpr.norm_list <- list(mmul_txi.gene.norm, mmul_txi.tx.norm, hg_txi.gene.norm, hg_txi.tx.norm)
datExpr.norm_name_list <- list('mmul_gene_post-normalizaton', 'mmul_tx_post-normalizaton', 'hg_gene_post-normalizaton', 'hg_tx_post-normalizaton')

# loop through each expression dataframe and determine which samples are beyond a connectivity z-score of 3
mapply(function(datExpr.norm, datExpr.norm_name) {
  
  # get connectivity z-scores for the normalized data
  normadj <- (0.5 + 0.5*bicor(datExpr.norm)^2)
  netsummary <- fundamentalNetworkConcepts(normadj)
  C <- netsummary$Connectivity
  Z.C <- (C-mean(C))/sqrt(var(C))
  
  # plot the connectivity z-score for the normalized data and set a cut-off at 3 STD
  plot(1:length(Z.C),Z.C,main=paste0("Outlier Plot of ", datExpr.norm_name),xlab = "\nSamples\n",ylab="\nConnectivity Z Score\n", ylim = c(-4, 1), col = datMeta$SeqBatch)
  legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = rev(as.numeric(as.factor(levels(datMeta$SeqBatch)))))
  abline(h = -3, col="red")
  
  # label the point that falls outside of the connectivity z-score
  text(which(names(Z.C) == '40482_V1'), Z.C[which(names(Z.C) == '40482_V1')], label = '40482_V1', pos=2, cex=0.8)
  
  # determine which samples fail the 'outlier threshold' via connectivity z-score
  print(which(Z.C < -3))
  
}, datExpr.norm_list, datExpr.norm_name_list)


# sample 40482_V1 failed the connectivity z-score test for both mmul_gene and hg_gene level analysis 
# so it was removed from all dataframes in the dataset
mmul_txi.gene <- mmul_txi.gene[,which(colnames(mmul_txi.gene) != '40482_V1')]
mmul_txi.tx <- mmul_txi.tx[,which(colnames(mmul_txi.tx) != '40482_V1')]
hg_txi.gene <- hg_txi.gene[,which(colnames(hg_txi.gene) != '40482_V1')]
hg_txi.tx <- hg_txi.tx[,which(colnames(hg_txi.tx) != '40482_V1')]
datMeta <- datMeta[which(rownames(datMeta) != '40482_V1'),]


### (5) Export dataframes to RData files #############################################################################

save(mmul_txi.gene, mmul_txi.tx, mmul_annotEns87, hg_txi.gene, hg_txi.tx, hg_annotEns84,
     datMeta, file = './working_data/filtered_unnormalized_datExpr_and_datMeta.RData')

dev.off()
dev.off()

rm(list=ls())
