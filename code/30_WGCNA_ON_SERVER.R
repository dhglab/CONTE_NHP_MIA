########################################################################################################################
#
#   30_WGCNA_ON_SERVER.R 
#
#   Perform WGCNA on 4 y/o non-human primate RNA-seq samples from DLPFC, ACC, HC, and V1. Conduct analysis via
#   multiple methods to determine the best way to perform WGCNA on this dataset.
#   
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

setwd('/u/home/n/npage98/project-geschwind/CONTE_NHP_MIA/WGCNA')

# Load libraries
library(WGCNA)
#enableWGCNAThreads(nThreads=8); # when requesting 8 CPUs
library(SummarizedExperiment)
library(edgeR)
library(SNFtool)
library(flashClust)
library(nlme)
library(gridExtra)
library(igraph) ## for plotting network plots

load('./filtered_unnormalized_datExpr_and_datMeta.RData')
load('./model_and_contrasts.RData')

#########################################################################################################################
#
#   (1) Prepare datExpr for WGCNA
#
#########################################################################################################################

source('./z_WGCNA_TOOLS.R')

datExpr <- mmul_txi.gene

# remove genes whose expression values are perfectly correlated/duplicates of eachother
datExpr <- datExpr[-which(duplicated(datExpr)),]


# get probe data for genes to be used for WGCNA from the annotation file
datProbes <- cbind(mmul_annotEns87$start_position[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$end_position[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$percentage_gc_content[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$ensembl_gene_id[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$hsapiens_homolog_associated_gene_name[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$hsapiens_homolog_ensembl_gene[match(rownames(datExpr), mmul_annotEns87$ensembl_gene_id)])
rownames(datProbes) <- rownames(datExpr)
colnames(datProbes) <- c('start_position', 'end_position', 'percentage_gc_content', 'ensembl_gene_id', 'external_gene_name', 'hsapiens_gene_id')
datProbes <- data.frame(datProbes)

datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr), method = 'TMM'), log = TRUE)

pdf('./30_WGCNA_ON_SERVER.pdf', 8, 8)

# visualize samples before regression
plot(cmdscale(dist(t(datExpr.norm))),col=datMeta$Region, main = 'Before Covariate Regression', xlab = 'MDS1', ylab = 'MDS2')
plot(cmdscale(dist(t(datExpr.norm))),col=datMeta$SeqBatch, main = 'Before Covariate Regression', xlab = 'MDS1', ylab = 'MDS2')

datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                      all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                      figureDir = './', figureName = 'mmul_gene_WGCNA')

# visualize samples after regression
plot(cmdscale(dist(t(datExpr.norm.reg))),col=datMeta$Region, main = 'After Covariate Regression', xlab = 'MDS1', ylab = 'MDS2')
plot(cmdscale(dist(t(datExpr.norm.reg))),col=datMeta$SeqBatch, main = 'After Covariate Regression', xlab = 'MDS1', ylab = 'MDS2')


#########################################################################################################################
#
#   (2) Perform WGCNA on all samples combined
#
#########################################################################################################################

sft = pickSoftThreshold(data = t(datExpr.norm.reg), networkType = "signed", corFnc = "bicor", verbose = 5, powerVector = seq(1, 30))
 
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n", main='', ylim = c(0, 1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = seq(1, 30), cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2",main='')
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = seq(1, 30), cex = 0.7, col="black")

dev.off()

powerEstimate <- 14

net <- blockwiseModules(datExpr = t(datExpr.norm.reg), power = powerEstimate, TOMType = 'signed', corType = 'bicor',
                        maxBlockSize = 15000, saveTOMFileBase = paste0('./WGCNA_sft_', powerEstimate, '_TOM'), saveTOMs = TRUE, verbose = 3)

load('./WGCNA_sft_14_TOM-block.1.RData')

save(net, file = paste0('./WGCNA_sft_', powerEstimate, '_network.RData'))
load(paste0('./WGCNA_sft_', powerEstimate, '_network.RData'))

geneTree = hclust(as.dist(1-TOM), method = "average")

# takes a long time to run
geneSigsColor <- relateDendroTraits(datExpr = datExpr.norm.reg, datMeta = datMeta, geneTree = geneTree, TOM = TOM, 
                                    outputName = paste0('WGCNA_sft_', powerEstimate, '_dendroTraits'))

finalDendroCut(datExpr = datExpr.norm.reg, geneTree = geneTree, dissTOM = TOM, geneSigsColor = geneSigsColor,
               softPower = 14, processedDataDir = paste0('./WGCNA_sft_', powerEstimate, '_Final'),
               mms = 100, ds = 4, dthresh = 0.1)

load(paste0('./WGCNA_sft_', powerEstimate, '_Final_Modules.RData'))

findEnrichment(MEs = MEs, datMeta = datMeta, datProbes = datProbes, moduleColors = moduleColors)











