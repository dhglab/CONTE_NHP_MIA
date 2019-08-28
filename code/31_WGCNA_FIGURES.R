########################################################################################################################
#
#   31_WGCNA_FIGURES.R 
#
#   Generate figures for WGCNA on 4 y/o non-human primate RNA-seq 
#   samples from DLPFC, ACC, HC, and V1
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(gridExtra)
library(ggplotify)
library(ggplot2)
library(ggrepel)
library(WGCNA)
library(edgeR)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')
load('./working_data/model_and_contrasts.RData')

load('./working_data/WGCNA/WGCNA_sft_14_dendroTraits.RData') 
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')
load('./working_data/WGCNA/WGCNA_sft_14_TOM-block.1.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Dendro tree for main figure
#   (2) Module trait association heatmap
#   (3) DE module cell-type enrichment
#   (4) Steelblue module boxplot
#   (5) Steelblue module hubs
#   (6) MEsteelblue g:Profiler enrichment
#   (7) PCA after regressing out top 4 seqPCs
#   (8) Show module correlation with region, gene length, and GC-content
#   (9) Show iteration through final dendroCut setting
#   (10) All module eigengene P-values
#
########################################################################################################################

### (1) Dendro tree for main figure ####################################################################################

mColorh <- t(rbind(mColorh1[,8], geneSigsColor))

colnames(mColorh) <- c('Modules', 'Combined MIA', '1st Trim MIA', '2nd Trim MIA')

mLabelh <- c('Modules', 'Combined MIA', '1st Trim MIA', '2nd Trim MIA')

pdf('./figures/Fig 2 Draft Pt 2.pdf', 2, 1.5)
plotDendroAndColors(geneTree, 
                    mColorh[,1], 
                    groupLabels = c(''),
                    addGuide = TRUE,
                    dendroLabels = FALSE,
                    ylab = NA,
                    main = NA,
                    cex.rowText = 0.6,
                    cex.colorLabels = 0.6,
                    cex.dendroLabels = 0.6,
                    marAll = c(0.2, 1, 0.2, 0.2))
dev.off()



### (2) Module trait association heatmap ###############################################################################

corfit <- duplicateCorrelation(t(MEs), mod, block = datMeta$Subject)
fit <- lmFit(t(MEs), mod, block = datMeta$Subject, correlation = corfit$consensus.correlation)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
results <- eBayes(fit_contrasts)

# get the results of the linear model
tt_all <- topTable(results, number = Inf, coef = 1, sort.by = 'none')
tt_DLPFC <- topTable(results, number = Inf, coef = 2, sort.by = 'none')
tt_ACC <- topTable(results, number = Inf, coef = 3, sort.by = 'none')
tt_HC <- topTable(results, number = Inf, coef = 4, sort.by = 'none')
tt_V1 <- topTable(results, number = Inf, coef = 5, sort.by = 'none')

tt_combined <- rbind(data.frame('module' = rownames(tt_DLPFC), 'logFC' = tt_DLPFC$logFC, 'p' = tt_DLPFC$P.Value, 'analysis' = '4'),
                     data.frame('module' = rownames(tt_ACC), 'logFC' = tt_ACC$logFC, 'p' = tt_ACC$P.Value, 'analysis' = '3'),
                     data.frame('module' = rownames(tt_HC), 'logFC' = tt_HC$logFC, 'p' = tt_HC$P.Value, 'analysis' = '2'),
                     data.frame('module' = rownames(tt_V1), 'logFC' = tt_V1$logFC, 'p' = tt_V1$P.Value, 'analysis' = '1'))

tt_combined$analysis <- factor(tt_combined$analysis, labels = rev(c('DLPFC', 'ACC', 'HC', 'V1')))

library(RColorBrewer)
palette <- colorRampPalette(rev(brewer.pal(11, 'RdBu')), space='Lab')

ME_trait_associations <- ggplot(tt_combined, aes(x = analysis, y = module, fill = -log10(p) * sign(logFC))) + 
  coord_flip() +
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(tt_combined$p < 0.05, '*', '')), 
            size = 8, color = 'yellow', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(-1.75, 1.75)) +
  ggtitle('Module Differential Expression') +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.margin=unit(c(2,5,5,2),"mm"),
        legend.box.margin=margin(-10,5,-10,-10))

save(ME_trait_associations, file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (3) DE module cell-type enrichment #################################################################################

library(EWCE)

load('./working_data/CellTypeData_PEC_Lake_Adult_Human_nucSeq.rda')

dim(ctd[[2]]$specificity)
head(ctd[[2]]$specificity)

load('./working_data/DGE_result_topTables.RData')
mmul_CON_vs_MIA_ALL <- DGE_result_topTables$mmul_CON_vs_MIA_ALL

set.seed(12211991)
reps = 10000
level1 = 2

gene2name <- match(names(moduleColors), mmul_annotEns87$ensembl_gene_id)
names(moduleColors) <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[gene2name]

MEsteelblue = names(moduleColors[which(moduleColors == 'steelblue')])
MEsalmon = names(moduleColors[which(moduleColors == 'salmon')])
MEpaleturquoise = names(moduleColors[which(moduleColors == 'paleturquoise')])

background_adult = rownames(ctd[[level1]]$specificity)
length(background_adult)
background_adult = background_adult[which(background_adult %in% mmul_CON_vs_MIA_ALL$external_gene_name)]
length(background_adult)

MEsteelblue = MEsteelblue[which(MEsteelblue %in% background_adult)]
MEsalmon = MEsalmon[which(MEsalmon %in% background_adult)]
MEpaleturquoise = MEpaleturquoise[which(MEpaleturquoise %in% background_adult)]

MEsteelblue_results = bootstrap.enrichment.test(sct_data=ctd, hits=MEsteelblue, bg=background_adult,reps=reps,annotLevel=level1,genelistSpecies="human",sctSpecies="human")
MEsalmon_results = bootstrap.enrichment.test(sct_data=ctd, hits=MEsalmon, bg=background_adult,reps=reps,annotLevel=level1,genelistSpecies="human",sctSpecies="human")
MEpaleturquoise_results = bootstrap.enrichment.test(sct_data=ctd, hits=MEpaleturquoise, bg=background_adult,reps=reps,annotLevel=level1,genelistSpecies="human",sctSpecies="human")

ewce_results <- rbind(data.frame(MEsteelblue_results$results[,c(1,3,5)], 'Module' = 'MEsteelblue'),
                      data.frame(MEsalmon_results$results[,c(1,3,5)], 'Module' = 'MEsalmon'),
                      data.frame(MEpaleturquoise_results$results[,c(1,3,5)], 'Module' = 'MEpaleturquoise'))

ewce_results$FDR <- p.adjust(ewce_results$p, method = 'fdr')
ewce_results$sd_from_mean <- ifelse(ewce_results$sd_from_mean > 0, ewce_results$sd_from_mean, 0)
print(ewce_results)

ewce_results$sig <- ifelse(ewce_results$FDR < 0.05, '*', '')

DE_module_cell_type_enrichment <- ggplot(data = ewce_results, aes(x = ewce_results$CellType, y = ewce_results$sd_from_mean, fill = ewce_results$Module)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_grid(cols = vars(Module)) +
  geom_text(aes(label=ewce_results$sig), vjust=0.5, color="red",
            position = position_dodge(0), size = 8) +
  scale_fill_manual(values = c('paleturquoise', 'salmon', 'steelblue')) +
  ggtitle('Cell-type Enrichment from PsychEncode') +
  ylab('Stdev from Mean') +
  ylim(c(0, 10)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(2,2,2,2),"mm"))

save(ME_trait_associations, DE_module_cell_type_enrichment,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (4) Steelblue module boxplot #######################################################################################

me_num <- which(colnames(MEs) == 'MEsteelblue')

MEsteelblue_boxplot <- ggplot(data.frame(MEs[,me_num]), aes(x=datMeta$MIA, y=MEs[,me_num])) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, size = 2, aes(color = datMeta$MIA)) +
  facet_grid(~ datMeta$Region) + 
  ggtitle(colnames(MEs)[me_num]) +
  theme_bw() +
  ylab(expression('Eigengene Expr')) +
  scale_color_discrete(labels = c("Control", "1st Trim MIA", "2nd Trim MIA")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 6),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(5, 2, 5, 2),"mm"))


# plot(MEsteelblue_boxplot)

save(ME_trait_associations, DE_module_cell_type_enrichment, MEsteelblue_boxplot,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (5) Steelblue module hubs ##########################################################################################

source('./code/z_WGCNA_TOOLS.R')

datExpr <- mmul_txi.gene
datExpr <- datExpr[-which(duplicated(datExpr)),]
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr), method = 'TMM'), log = TRUE)
datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                           all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                           figureDir = './docs/trait_correlations/', figureName = 'mmul_gene_WGCNA')


library(igraph)

kME = signedKME(t(datExpr.norm.reg), MEs, corFnc = "bicor")

hubGenes <- rownames(kME)[order(kME[,which(colnames(kME) == 'kMEsteelblue')], decreasing = T)[1:20]]
hubGene.symbols <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[match(hubGenes, mmul_annotEns87$ensembl_gene_id)]

hubGene.symbols = ifelse(mmul_annotEns87$hsapiens_homolog_associated_gene_name[match(hubGenes, mmul_annotEns87$ensembl_gene_id)] != '',
                         mmul_annotEns87$hsapiens_homolog_associated_gene_name[match(hubGenes, mmul_annotEns87$ensembl_gene_id)],
                         mmul_annotEns87$ensembl_gene_id[match(hubGenes, mmul_annotEns87$ensembl_gene_id)])

adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "bicor", power=10)
adjMat[adjMat < quantile(adjMat,0.1)]=0
graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
plotcord= data.frame(layout_with_dh(graph))
colnames(plotcord) = c("X1","X2")
edgelist <- get.edgelist(graph,names = F)
edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
colnames(edges) <- c("X1","Y1","X2","Y2")
plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))


MEsteelblue_hubs <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edges, size = 0.5, color = "grey") + 
  geom_point(aes(X1, X2), data = plotcord, color='steelblue', size = 3) + 
  geom_text(aes(x = X1, y = X2 + 0.2, label = gene, fontface = 'bold'), 
                  size = 1.8, data = plotcord, show.legend = FALSE, vjust = 'inward',hjust = 'inward') +
  theme_classic() + labs(x="", y="") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size = 6),
        plot.margin=unit(c(0, 2, 0, 0), "mm"))

plot(MEsteelblue_hubs)

save(ME_trait_associations, DE_module_cell_type_enrichment, MEsteelblue_boxplot,
     MEsteelblue_hubs,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (6) MEsteelblue g:Profiler enrichment ##############################################################################

library(gProfileR)

# important to reload the moduleColors so that we can get gene IDs
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')

query_background <- rownames(mmul_CON_vs_MIA_ALL)
#query_background <- query_background[which(query_background != '')]

MEsteelblue_ids <- names(moduleColors)[which(moduleColors == 'steelblue')]
MEsteelblue_kMEs <- kME$kMEsteelblue[which(rownames(kME) %in% MEsteelblue_ids)]
names(MEsteelblue_kMEs) <- rownames(kME)[which(rownames(kME) %in% MEsteelblue_ids)]
MEsteelblue_kMEs <- MEsteelblue_kMEs[order(MEsteelblue_kMEs, decreasing = TRUE)]

go = gprofiler(names(MEsteelblue_kMEs), organism="mmulatta", custom_bg = query_background,
               correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
               region_query = F,max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
               include_graph = F,src_filter = c("GO", "KEGG"))

go <- go[order(go$p.value, decreasing = FALSE),]
GO_summary <- data.frame('term' = go$term.name, 'p' = go$p.value, 'core_enrichment' = go$intersection, stringsAsFactors = FALSE)
GO_summary <- GO_summary[c(1:5),]
GO_summary$term <- factor(GO_summary$term,
                          levels = GO_summary$term,
                          labels = c('steroid metabolic process',
                                     'peptidase inhibitor activity',
                                     'extracellular space',
                                     'maintenance of mitochondrion\nlocation',
                                     'SREBP-SCAP-Insig complex'))

MEsteelblue_GO <- ggplot(data = GO_summary, aes(x = GO_summary$term, y = -log10(GO_summary$p))) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_summary$term)) +
  #facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  ylim(c(0, 1.8)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 10, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 5),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(-2,5,5,-1),"mm"))

plot(MEsteelblue_GO)

save(ME_trait_associations, DE_module_cell_type_enrichment, MEsteelblue_boxplot,
     MEsteelblue_hubs, MEsteelblue_GO,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (7) PCA after regressing out top 4 seqPCs ##########################################################################

ePCs <- prcomp(t(scale(t(datExpr.norm.reg), scale = TRUE)), center = FALSE)
ePCs <- data.frame(ePCs$rotation)

ePC_reg_plot <- ggplot(ePCs, aes(x = PC1, y = PC2, color = datMeta$Region, shape = datMeta$SeqBatch, fill = datMeta$MIA)) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(22, 21), labels = c('Batch 1', 'Batch 2')) +
  scale_fill_manual(values = c('white', 'grey'), labels = c('Batch 1', 'Batch 2')) +
  ggtitle('Post SeqPC Regression') +
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

plot(ePC_reg_plot)

save(ME_trait_associations, DE_module_cell_type_enrichment, MEsteelblue_boxplot,
     MEsteelblue_hubs, MEsteelblue_GO, ePC_reg_plot,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')



### (8) Show module correlation with region, gene length, and GC-content ###############################################

library(scales)

# get gene lengths and GC-content
datProbes <- cbind(mmul_annotEns87$end_position[match(rownames(datExpr.norm.reg), mmul_annotEns87$ensembl_gene_id)] - mmul_annotEns87$start_position[match(rownames(datExpr.norm.reg), mmul_annotEns87$ensembl_gene_id)],
                   mmul_annotEns87$percentage_gc_content[match(rownames(datExpr.norm.reg), mmul_annotEns87$ensembl_gene_id)])
rownames(datProbes) <- rownames(datExpr.norm.reg)
colnames(datProbes) <- c('gene_length', 'percentage_gc_content')
datProbes <- data.frame(datProbes)


geneSigs=matrix(NA,nrow = 4,ncol=nrow(datExpr.norm.reg)) # create a vector to hold the data

for( i in 1:ncol(geneSigs)) {
  # calculate r correlation value for numeric variables
  # calculate adjusted R^2s square-root for categorical variables (factors)
  exprvec=as.numeric(datExpr.norm.reg[i,]) # get the expression vector for ith gene
  
  DLPFC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="DLPFC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="DLPFC"))))$coefficients[2,1]) # for DLPFC
  ACC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="ACC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="ACC"))))$coefficients[2,1]) # for ACC
  HC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="HC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="HC"))))$coefficients[2,1]) # for HC
  V1_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="V1"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="V1"))))$coefficients[2,1]) # for V1
  
  geneSigs[,i]=c(DLPFC_r, ACC_r, HC_r, V1_r)
}

colnames(geneSigs) = rownames(datExpr.norm.reg)

# center and scale gene length and GC-content
datProbes$gene_length <- rescale(log2(datProbes$gene_length), to = c(0, 1))
datProbes$percentage_gc_content <- rescale(datProbes$percentage_gc_content, to = c(0, 1))

geneSigs <- rbind(geneSigs, t(datProbes$gene_length), t(datProbes$percentage_gc_content))


geneSigsColor=matrix(NA,nrow=nrow(geneSigs),ncol=nrow(datExpr)) # create a vector to hold the data
for ( i in 1:nrow(geneSigsColor)) {
  geneSigsColor[i,] = numbers2colors(as.numeric(geneSigs[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 
}




rownames(geneSigsColor) = c("DLPFC", "ACC", "HC", "V1", "Length", '%GC')
colnames(geneSigsColor) = colnames(geneSigs)

mColorh2=cbind(mColorh[,1],t(geneSigsColor))
mLabelh2=c(mLabelh[1],rownames(geneSigsColor))

pdf('./figures/Fig S4 Draft Pt 3.pdf', 4.5, 2.5)
plotDendroAndColors(geneTree, 
                    mColorh2, 
                    groupLabels = mLabelh2,
                    addGuide = TRUE,
                    dendroLabels = FALSE,
                    ylab = 'Height',
                    main = 'WGCNA',
                    cex.rowText = 0.6,
                    cex.colorLabels = 0.6,
                    cex.dendroLabels = 0.6,
                    marAll = c(0.5, 4, 1, 0))
dev.off()



### (9) Show iteration through final dendroCut setting #################################################################

library(stringr)
mLabelh1 <- mLabelh1[1:18]
mLabelh1 <- str_remove(mLabelh1, '\n ')
mLabelh1 <- str_replace(mLabelh1, '= ', '=')
mLabelh1 <- str_replace(mLabelh1, '= 0', '=0')
mLabelh1 <- str_replace(mLabelh1, '  ', ';  ')
mLabelh1 <- str_replace(mLabelh1, '  d', ';  d')

pdf('./figures/Fig S4 Draft Pt 4.pdf', 3, 3.5)
plotDendroAndColors(geneTree, 
                    mColorh1[,1:18], 
                    groupLabels = mLabelh1[1:18],
                    addGuide = TRUE,
                    dendroLabels = FALSE,
                    cex.rowText = 0.2,
                    cex.colorLabels = 0.3,
                    cex.dendroLabels = 0.2,
                    marAll = c(0.5, 4, 1, 0))
dev.off()



### (10) All module eigengene P-values #################################################################################

corfit <- duplicateCorrelation(t(MEs), mod, block = datMeta$Subject)
fit <- lmFit(t(MEs), mod, block = datMeta$Subject, correlation = corfit$consensus.correlation)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
results <- eBayes(fit_contrasts)

# get the results of the linear model
tt_all <- topTable(results, number = Inf, coef = 1, sort.by = 'none')
tt_DLPFC <- topTable(results, number = Inf, coef = 2, sort.by = 'none')
tt_ACC <- topTable(results, number = Inf, coef = 3, sort.by = 'none')
tt_HC <- topTable(results, number = Inf, coef = 4, sort.by = 'none')
tt_V1 <- topTable(results, number = Inf, coef = 5, sort.by = 'none')

tt_combined <- rbind(data.frame('module' = rownames(tt_all), 'logFC' = tt_all$logFC, 'p' = tt_all$P.Value, 'analysis' = '1'),
                     data.frame('module' = rownames(tt_DLPFC), 'logFC' = tt_DLPFC$logFC, 'p' = tt_DLPFC$P.Value, 'analysis' = '2'),
                     data.frame('module' = rownames(tt_ACC), 'logFC' = tt_ACC$logFC, 'p' = tt_ACC$P.Value, 'analysis' = '3'),
                     data.frame('module' = rownames(tt_HC), 'logFC' = tt_HC$logFC, 'p' = tt_HC$P.Value, 'analysis' = '4'),
                     data.frame('module' = rownames(tt_V1), 'logFC' = tt_V1$logFC, 'p' = tt_V1$P.Value, 'analysis' = '5'))

tt_combined$analysis <- factor(tt_combined$analysis, labels = c('Combined', 'DLPFC', 'ACC', 'HC', 'V1'))

library(RColorBrewer)
palette <- colorRampPalette(rev(brewer.pal(11, 'RdBu')), space='Lab')

All_ME_trait_associations <- ggplot(tt_combined, aes(x = analysis, y = module, fill = -log10(p) * sign(logFC))) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = paste0('p =\n', round(tt_combined$p, 2))), 
            size = 1.8, color = ifelse(tt_combined$p < 0.2, 'white', 'black'), fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(-1.75, 1.75)) +
  ggtitle('Module Differential Expression') +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-10,0,0,0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.margin=unit(c(2,5,10,2),"mm"),
        legend.box.margin=margin(-10,5,-10,-10))

plot(All_ME_trait_associations)

save(ME_trait_associations, DE_module_cell_type_enrichment, MEsteelblue_boxplot,
     MEsteelblue_hubs, MEsteelblue_GO, ePC_reg_plot, All_ME_trait_associations,
     file = './working_data/plots/31_WGCNA_FIGURES.RData')













