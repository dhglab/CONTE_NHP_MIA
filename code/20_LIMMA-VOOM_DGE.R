########################################################################################################################
#
#   20_LIMMA-VOOM_DGE.R 
#
#   Perform differential gene expression on 4 y/o non-human primate RNA-seq samples from DLPFC, ACC, HC, and V1
#   Conduct analysis on all samples combined and on samples split by brain region
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(gridExtra)
library(ggplot2)
library(ggrepel)
library(limma)
library(edgeR)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')

# create dataframe to load topTable results into
DGE_result_topTables <- list()

########################################################################################################################
#
#   Pipeline:
#
#   (1) Calculate DGE on all samples combined by CON vs. MIA
#   (2) Calculate DGE on all samples combined by Group
#   (3) Save topTables as a list for downstream analysis
#
########################################################################################################################

### (1) Calculate DGE on all samples combined by CON vs. MIA ###########################################################

# create model matrix for DGE and include MG_seqPCs 1-4
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

datExpr_list <- list(mmul_txi.gene, hg_txi.gene)
datExpr_name_list <- list('mmul', 'hg')
analysis <- list('CON_vs_MIA_ALL', 'CON_vs_MIA_DLPFC', 'CON_vs_MIA_ACC', 'CON_vs_MIA_HC', 'CON_vs_MIA_V1')

pdf('./docs/20_CONTE_NHP_MIA_DGE_pt1.pdf', 8, 8)

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
    
    # add external gene names!!
    if(datExpr_name == 'mmul') {
      
      gene2name <- match(rownames(tt.fullModel), mmul_annotEns87$ensembl_gene_id)
      tt.fullModel$external_gene_name <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[gene2name]
      
    } else {
      
      gene2name <- match(rownames(tt.fullModel), hg_annotEns84$ensembl_gene_id)
      tt.fullModel$external_gene_name <- hg_annotEns84$external_gene_name[gene2name]
      
    }
    
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
      geom_text_repel(data=tt.fullModel[tt.fullModel$P.Val < nomPVALcutoff,], aes(label=external_gene_name),col='grey60') +   #Add text labels for 'significant' genes only
      theme_bw() + ggtitle(paste0(datExpr_name, " ", analysis[coef_num], " DGE Volcano Plot (after fdrtool"))
    
    plot(p2)
    
    # save topTable for each analysis to a dataframe that can be referenced later (note: need to use '<<-' so assignment stays outside of the loop)
    DGE_result_topTables[[paste0(datExpr_name, '_', analysis[[coef_num]])]] <<- tt.fullModel
    
  }
  
}, datExpr_list, datExpr_name_list)

dev.off()



### (2) Calculate DGE on all samples combined by Group #################################################################

mod <- cbind(ifelse(datMeta$Group == 'con' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$Group == 'con' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$Group == 'con' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$Group == 'con' & datMeta$Region == 'V1', 1, 0),
             ifelse(datMeta$Group == 'poly1' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$Group == 'poly1' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$Group == 'poly1' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$Group == 'poly1' & datMeta$Region == 'V1', 1, 0),
             ifelse(datMeta$Group == 'poly2' & datMeta$Region == 'DLPFC', 1, 0),
             ifelse(datMeta$Group == 'poly2' & datMeta$Region == 'ACC', 1, 0),
             ifelse(datMeta$Group == 'poly2' & datMeta$Region == 'HC', 1, 0),
             ifelse(datMeta$Group == 'poly2' & datMeta$Region == 'V1', 1, 0),
             datMeta$MG_seqPC1,
             datMeta$MG_seqPC2,
             datMeta$MG_seqPC3,
             datMeta$MG_seqPC4)

colnames(mod) <- c('con_DLPFC', 'con_ACC', 'con_HC', 'con_V1',
                   'poly1_DLPFC', 'poly1_ACC', 'poly1_HC', 'poly1_V1',
                   'poly2_DLPFC', 'poly2_ACC', 'poly2_HC', 'poly2_V1',
                   'MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4')

rownames(mod) <- rownames(datMeta)

contrast.matrix <- makeContrasts((1/4) * (poly1_DLPFC + poly1_ACC + poly1_HC + poly1_V1 -
                                            con_DLPFC - con_ACC - con_HC - con_V1),
                                 (1/4) * (poly2_DLPFC + poly2_ACC + poly2_HC + poly2_V1 -
                                            con_DLPFC - con_ACC - con_HC - con_V1),
                                 (1/4) * (poly2_DLPFC + poly2_ACC + poly2_HC + poly2_V1 -
                                            poly1_DLPFC - poly1_ACC - poly1_HC - poly1_V1),
                                 poly1_DLPFC - con_DLPFC,
                                 poly2_DLPFC - con_DLPFC,
                                 poly2_DLPFC - poly1_DLPFC,
                                 poly1_ACC - con_ACC,
                                 poly2_ACC - con_ACC,
                                 poly2_ACC - poly1_ACC,
                                 poly1_HC - con_HC,
                                 poly2_HC - con_HC,
                                 poly2_HC - poly1_HC,
                                 poly1_V1 - con_V1,
                                 poly2_V1 - con_V1,
                                 poly2_V1 - poly1_V1,
                                 levels = mod)

head(mod)


datExpr_list <- list(mmul_txi.gene, hg_txi.gene)
datExpr_name_list <- list('mmul', 'hg')
analysis <- list('con_vs_poly1', 'con_vs_poly2', 'poly1_vs_poly2',
                 'con_vs_poly1_DLPFC', 'con_vs_poly2_DLPFC', 'poly1_vs_poly2_DLPFC',
                 'con_vs_poly1_ACC', 'con_vs_poly2_ACC', 'poly1_vs_poly2_ACC',
                 'con_vs_poly1_HC', 'con_vs_poly2_HC', 'poly1_vs_poly2_HC',
                 'con_vs_poly1_V1', 'con_vs_poly2_V1', 'poly1_vs_poly2_V1')



pdf('./docs/20_CONTE_NHP_MIA_DGE_pt2.pdf', 8, 8)

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
  for(coef_num in seq(1:15)){
    
    # Top table = differential expression summary statistics for each gene
    tt.fullModel = topTable(lm,coef=coef_num, number=Inf, sort.by = 'P')
    tt.fullModel$FDR = fdrtool::fdrtool(tt.fullModel$t, plot = FALSE)$qval
    
    # add external gene names!!
    if(datExpr_name == 'mmul') {
      
      gene2name <- match(rownames(tt.fullModel), mmul_annotEns87$ensembl_gene_id)
      tt.fullModel$external_gene_name <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[gene2name]
      
    } else {
      
      gene2name <- match(rownames(tt.fullModel), hg_annotEns84$ensembl_gene_id)
      tt.fullModel$external_gene_name <- hg_annotEns84$external_gene_name[gene2name]
      
    }
    
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
      geom_text_repel(data=tt.fullModel[tt.fullModel$P.Val < nomPVALcutoff,], aes(label=external_gene_name),col='grey60') +   #Add text labels for 'significant' genes only
      theme_bw() + ggtitle(paste0(datExpr_name, " ", analysis[coef_num], " DGE Volcano Plot (after fdrtool"))
    
    plot(p2)
    
    # save topTable for each analysis to a dataframe that can be referenced later (note: need to use '<<-' so assignment stays outside of the loop)
    DGE_result_topTables[[paste0(datExpr_name, '_', analysis[[coef_num]])]] <<- tt.fullModel
    
  }
  
}, datExpr_list, datExpr_name_list)

dev.off()



### (3) Save topTables as a list for downstream analysis ###############################################################

save(DGE_result_topTables, file = './working_data/DGE_result_topTables.RData')



