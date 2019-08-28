########################################################################################################################
#
#   40_TRANSPOSABLE_ELEMENTS_DE.R 
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

library(ggplot2)
library(ggrepel)
library(edgeR)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Gather TE data from RepEnrich and perfor QC
#   (2) Perform global TE differential expression with Limma-voom
#   (3) Perform class and family transposable element DE with a linear model
#
########################################################################################################################

### (1) Gather TE data from RepEnrich and perfor QC ####################################################################

pdf('./docs/40_TRANSPOSABLE_ELEMENTS_DE_pt1.pdf', 8, 8)

rep_datExpr <- read.csv(file = './raw_data/TransposableElements/RepEnrich2/TE_fraction_counts.csv', header = TRUE, sep = '\t')
rownames(rep_datExpr) <- rep_datExpr$TE
rep_datExpr <- rep_datExpr[,c(2:53)]
colnames(rep_datExpr) <- substring(colnames(rep_datExpr), 2)


class_rep_datExpr <- read.csv(file = './raw_data/TransposableElements/RepEnrich2/TE_class_fraction_counts.csv', header = TRUE, sep = '\t')
rownames(class_rep_datExpr) <- class_rep_datExpr$TE
class_rep_datExpr <- class_rep_datExpr[,c(2:53)]
colnames(class_rep_datExpr) <- substring(colnames(class_rep_datExpr), 2)


family_rep_datExpr <- read.csv(file = './raw_data/TransposableElements/RepEnrich2/TE_family_fraction_counts.csv', header = TRUE, sep = '\t')
rownames(family_rep_datExpr) <- family_rep_datExpr$TE
family_rep_datExpr <- family_rep_datExpr[,c(2:53)]
colnames(family_rep_datExpr) <- substring(colnames(family_rep_datExpr), 2)

# remove outlier V1 sample from all dataframes
rep_datExpr <- rep_datExpr[,which(colnames(rep_datExpr) != '40482_V1')]
class_rep_datExpr <- class_rep_datExpr[,which(colnames(class_rep_datExpr) != '40482_V1')]
family_rep_datExpr <- family_rep_datExpr[,which(colnames(family_rep_datExpr) != '40482_V1')]


# before filtering low expressed TEs
rep_datExpr.norm <- cpm(calcNormFactors(DGEList(rep_datExpr), method = 'TMM'), log = TRUE)
i=1; plot(density((rep_datExpr.norm[,i]),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
          main = 'Transposable Elements (Pre-filter)', xlab="cpm((TMM normalized counts), log = TRUE)") #, xlim=c(-5,20), ylim=c(0,0.20))
for(i in 2:dim(rep_datExpr.norm)[2]){
  lines(density((rep_datExpr.norm[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
}
legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = as.numeric(as.factor(levels(datMeta$SeqBatch))))

# remove low expresssed TEs
genes_to_keep <- filterByExpr(DGEList(rep_datExpr)) #need to mention this in the methods section
table(genes_to_keep)
rep_datExpr <- rep_datExpr[genes_to_keep,] # removed 72 low expressed TEs

# after filtering low expressed TEs
rep_datExpr.norm <- cpm(calcNormFactors(DGEList(rep_datExpr), method = 'TMM'), log = TRUE)
i=1; plot(density((rep_datExpr.norm[,i]),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
          main = 'Transposable Elements (Post-filter)', xlab="cpm((TMM normalized counts), log = TRUE)") #, xlim=c(-5,20), ylim=c(0,0.20))
for(i in 2:dim(rep_datExpr.norm)[2]){
  lines(density((rep_datExpr.norm[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
}
legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = as.numeric(as.factor(levels(datMeta$SeqBatch))))

# after filtering low expressed TEs and TMM normalization
rep_datExpr.norm <- cpm(calcNormFactors(DGEList(rep_datExpr), method = 'TMM'), log = TRUE)
i=1; plot(density((rep_datExpr.norm[,i]),na.rm=T),col = as.numeric(datMeta$SeqBatch)[i],
          main = 'Transposable Elements (Post-filter)', xlab="cpm((TMM normalized counts), log = TRUE)") #, xlim=c(-5,20), ylim=c(0,0.20))
for(i in 2:dim(rep_datExpr.norm)[2]){
  lines(density((rep_datExpr.norm[,i]),na.rm=T), col = as.numeric(datMeta$SeqBatch)[i],)
}
legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = as.numeric(as.factor(levels(datMeta$SeqBatch))))

# plot first two expression PCs
TE_ePCs <- prcomp(t(scale(t(rep_datExpr.norm), scale = TRUE)), center = FALSE)
TE_ePCs <- data.frame(TE_ePCs$rotation)
TE_ePC_plot <- ggplot(TE_ePCs, aes(x = PC1, y = PC2, color = datMeta$Region, shape = datMeta$SeqBatch, fill = datMeta$MIA)) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(22, 21), labels = c('Batch 1', 'Batch 2')) +
  scale_fill_manual(values = c('white', 'grey'), labels = c('Batch 1', 'Batch 2')) +
  ggtitle('Quality Control') +
  labs(color = 'Region', shape = 'Ctrl', fill = 'MIA') + 
  theme_bw() +
  guides(color = guide_legend(order = 0),
         shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes = list(shape = c(22, 21), fill = c('grey', 'grey')))) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(2.5, 'mm'),
        legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.direction = 'vertical',
        legend.box.margin = margin(-10, -5, -10, -10))

plot(TE_ePC_plot)



### (2) Perform global TE differential expression with Limma-voom ######################################################

mod <- cbind(ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'V1', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'V1', 1, 0),
             datMeta$MG_seqPC1,
             datMeta$MG_seqPC2,
             datMeta$MG_seqPC3,
             datMeta$MG_seqPC4)

# rename columns of the model matrix
colnames(mod) <- c('CON_DLPFC', 'CON_ACC', 'CON_HC', 'CON_V1',
                   'MIA_DLPFC', 'MIA_ACC', 'MIA_HC', 'MIA_V1',
                   'MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4')

# copy the row names for the model matrix from the metadata
rownames(mod) <- rownames(datMeta)

# create a contrast table for each test of differential gene expression we interested in
contrast.matrix <- makeContrasts((1/4) * (MIA_DLPFC + MIA_ACC + MIA_HC + MIA_V1 -
                                            CON_DLPFC - CON_ACC - CON_HC - CON_V1),
                                 MIA_DLPFC - CON_DLPFC,
                                 MIA_ACC - CON_ACC,
                                 MIA_HC - CON_HC,
                                 MIA_V1 - CON_V1, levels = mod)

head(mod)


datExpr_list <- list(rep_datExpr)
datExpr_name_list <- list('RepEnrich')
analysis <- list('CON_vs_MIA_ALL', 'CON_vs_MIA_DLPFC', 'CON_vs_MIA_ACC', 'CON_vs_MIA_HC', 'CON_vs_MIA_V1')

TE_result_topTables <- list()

mapply(function(datExpr, datExpr_name) {
  
  datExpr.counts <- datExpr
  dim(datExpr.counts)
  
  counts = DGEList(datExpr.counts) # Recomputing with outilers removed
  counts = calcNormFactors(counts,method = 'TMM')  #Perform Trimmed Mean of M-value (TMM) library size normalization 
  voom <- voom(design = mod, counts, plot = TRUE) # VOOM normalization for Limma, see: https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
  corfit= duplicateCorrelation(voom,mod,block=datMeta$Subject)   # This is for random effect repeated measure; see: https://support.bioconductor.org/p/59700/
  voom = voom(design= mod, counts, block=datMeta$Subject, correlation = corfit$consensus, plot=T) # Repeat with random effect blocking
  corfit= duplicateCorrelation(voom,mod,block=datMeta$Subject)   # This is for random effect repeated measure; see: https://support.bioconductor.org/p/59700/
  
  lm = lmFit(object = voom, design = mod, block=datMeta$Subject, correlation = corfit$consensus)   #Fits the linear model
  lm <- contrasts.fit(lm, contrast.matrix)
  lm = eBayes(lm)   # Empiric bayes moderated T-statistics
  
  # loops through each coefficent in the contrast matrix to generate a volcano plot
  for(coef_num in seq(1:5)){
    
    # Top table = differential expression summary statistics for each gene
    tt.fullModel = topTable(lm,coef=coef_num, number=Inf, sort.by = 'P')
    tt.fullModel$FDR = fdrtool::fdrtool(tt.fullModel$t, plot = FALSE)$qval
    
    tt.fullModel['external_gene_name'] <- rownames(tt.fullModel)
    
    hist(tt.fullModel$P.Value) # P value distribution
    table(tt.fullModel$adj.P.Val<.1)
    
    # Volcano plot of effect size vs -log(P value)
    nomPVALcutoff=0.005
    p1 <- ggplot(tt.fullModel, aes(x=logFC, y=-log10(P.Value), color = adj.P.Val < .1)) + geom_point() + #xlim(-2,2) +
      geom_text_repel(data=tt.fullModel[tt.fullModel$P.Val < nomPVALcutoff,], aes(label=external_gene_name),col='grey60') +   #Add text labels for 'significant' genes only
      theme_bw() + ggtitle(paste0(datExpr_name, " ", analysis[coef_num], " DGE Volcano Plot"))
    
    plot(p1)
    
    # Volcano plot of effect size vs -log(P value)
    nomPVALcutoff=0.005
    p2 <- ggplot(tt.fullModel, aes(x=logFC, y=-log10(P.Value), color = FDR < .1)) + geom_point() + #xlim(-2,2) +
      geom_text_repel(data=tt.fullModel[tt.fullModel$P.Val < nomPVALcutoff | tt.fullModel$FDR < 0.1,], aes(label=external_gene_name),col='grey60') +   #Add text labels for 'significant' genes only
      theme_bw() + ggtitle(paste0(datExpr_name, " ", analysis[coef_num], " DGE Volcano Plot (after fdrtool"))
    
    plot(p2)
    
    # save topTable for each analysis to a dataframe that can be referenced later (note: need to use '<<-' so assignment stays outside of the loop)
    TE_result_topTables[[paste0(datExpr_name, '_', analysis[[coef_num]])]] <<- tt.fullModel
    
  }
  
}, datExpr_list, datExpr_name_list)

dev.off()



### (3) Perform class and family transposable element DE with a linear model ###########################################

family_rep_datExpr.norm <- cpm(calcNormFactors(DGEList(family_rep_datExpr), method = 'TMM'), log = TRUE)
dim(family_rep_datExpr)

mod <- cbind(ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$MIA == 'CON' & datMeta$Region == 'V1', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$MIA == 'MIA' & datMeta$Region == 'V1', 1, 0),
             datMeta$MG_seqPC1,
             datMeta$MG_seqPC2,
             datMeta$MG_seqPC3,
             datMeta$MG_seqPC4)

# rename columns of the model matrix
colnames(mod) <- c('CON_DLPFC', 'CON_ACC', 'CON_HC', 'CON_V1',
                   'MIA_DLPFC', 'MIA_ACC', 'MIA_HC', 'MIA_V1',
                   'MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4')

# copy the row names for the model matrix from the metadata
rownames(mod) <- rownames(datMeta)

# create a contrast table for each test of differential gene expression we interested in
contrast.matrix <- makeContrasts((1/4) * (MIA_DLPFC + MIA_ACC + MIA_HC + MIA_V1 -
                                            CON_DLPFC - CON_ACC - CON_HC - CON_V1),
                                 MIA_DLPFC - CON_DLPFC,
                                 MIA_ACC - CON_ACC,
                                 MIA_HC - CON_HC,
                                 MIA_V1 - CON_V1, levels = mod)



# fits a linear model to the data that will be used to determine module enrichment
fit <- lmFit(family_rep_datExpr, mod)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
results <- eBayes(fit_contrasts)


# get the results of the linear model
tt_all <- topTable(results, number = Inf, coef = 1, sort.by = 'none')
tt_DLPFC <- topTable(results, number = Inf, coef = 2, sort.by = 'none')
tt_ACC <- topTable(results, number = Inf, coef = 3, sort.by = 'none')
tt_HC <- topTable(results, number = Inf, coef = 4, sort.by = 'none')
tt_V1 <- topTable(results, number = Inf, coef = 5, sort.by = 'none')

datExpr.norm <- family_rep_datExpr.norm

pdf('./docs/40_TRANSPOSABLE_ELEMENTS_DE_pt2.pdf', 8, 6)
for(i in 1:dim(family_rep_datExpr)[1]){
  
  p_list <- paste0('Combined: ', round(tt_all$P.Value[i], 3), '\n',
                   'DLPFC: ', round(tt_DLPFC$P.Value[i], 3), '\n',
                   'ACC: ', round(tt_ACC$P.Value[i], 3), '\n',
                   'HC: ', round(tt_HC$P.Value[i], 3), '\n',
                   'V1: ', round(tt_V1$P.Value[i], 3), '\n')
  
  gene_id <- rownames(family_rep_datExpr.norm)[i]
  
  p1 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_id),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_id),], color = datMeta$Group)) + 
    guides(color=guide_legend(title="Group")) + 
    geom_boxplot(outlier.shape = NA, lwd = 1.2) +
    #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
    geom_point(shape = 16, size = 4, position = position_dodge(width = 0.75)) +
    facet_grid(. ~ datMeta$Region) + 
    ggtitle(paste0(gene_id, '\n', p_list)) +
    theme_bw() +
    ylab(expression('Expr (log '[2]*' (CPM))')) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_text(size = 18, angle = 45, margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x=element_blank(),
          strip.text.x = element_text(size = 20),
          legend.title=element_text(size=12),
          legend.text=element_text(size=16))
  
  plot(p1)
  
}
dev.off()

save(TE_result_topTables, file = './working_data/TE_result_topTables.RData')



















