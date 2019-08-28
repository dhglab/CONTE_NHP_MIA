########################################################################################################################
#
#   32_WGCNA_FIGURES.R 
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
load('./working_data/DGE_result_topTables.RData')
load('./working_data/model_and_contrasts.RData')

mmul_CON_vs_MIA_ALL <- DGE_result_topTables$mmul_CON_vs_MIA_ALL

load('./working_data/WGCNA/WGCNA_sft_14_dendroTraits.RData') 
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')
#load('./working_data/WGCNA/WGCNA_sft_14_TOM-block.1.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Paleturquoise module boxplot
#   (2) Paleturquoise module hubs
#   (3) MEPaleturquoise g:Profiler enrichment
#   (4) Salmon module boxplot
#   (5) Salmon module hubs
#   (6) MEsalmon g:Profiler enrichment
#   (7) Heatmaps of cell-type specificity for MEpaleturquoise
#   (8) Turquoise module boxplot
#   (9) Turquoise module hubs
#   (10) MEturquoise g:Profiler enrichment
#   (11) MEturquoise cell-type enrichment
#   (12) Combined WGCNA module DE plot
#
########################################################################################################################

MEsalmon = names(moduleColors[which(moduleColors == 'salmon')])
MEpaleturquoise = names(moduleColors[which(moduleColors == 'paleturquoise')])
MEturquoise = names(moduleColors[which(moduleColors == 'turquoise')])

HC_MEs <- MEs[which(datMeta$Region == 'HC'),]
HC_datMeta <- datMeta[which(datMeta$Region == 'HC'),]

### (1) Paleturquoise module boxplot #######################################################################################

i <- which(colnames(MEs) == 'MEpaleturquoise')

MEpaleturquoise_boxplot <- ggplot(data.frame(HC_MEs[,i]), aes(x=HC_datMeta$MIA, y=HC_MEs[,i])) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, size = 2, aes(color = HC_datMeta$MIA)) +
  facet_grid(~ HC_datMeta$Region) + 
  ggtitle(colnames(HC_MEs)[i]) +
  theme_bw() +
  ylab(expression('Eigengene Expr')) +
  scale_color_discrete(labels = c("Control", "1st Trim MIA", "2nd Trim MIA")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(2,1,2,1),"mm"))


# plot(MEpaleturquoise_boxplot)

save(MEpaleturquoise_boxplot, file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (2) Paleturquoise module hubs ##########################################################################################

source('./code/z_WGCNA_TOOLS.R')

datExpr <- mmul_txi.gene
datExpr <- datExpr[-which(duplicated(datExpr)),]
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr), method = 'TMM'), log = TRUE)
datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                           all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                           figureDir = './docs/trait_correlations/', figureName = 'mmul_gene_WGCNA')


library(igraph)

kME = signedKME(t(datExpr.norm.reg), MEs, corFnc = "bicor")

hubGenes <- rownames(kME)[order(kME[,which(colnames(kME) == 'kMEpaleturquoise')], decreasing = T)[1:20]]
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
MEpaleturquoise_plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))


MEpaleturquoise_hubs <- ggplot(data = MEpaleturquoise_plotcord) + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edges, size = 0.5, color = "grey") + 
  geom_point(aes(X1, X2), data = MEpaleturquoise_plotcord, color='paleturquoise', size = 3) + 
  geom_text_repel(aes(x = X1, y = X2 + 0.1, label = gene, fontface = 'bold'), 
                  size = 2, data = MEpaleturquoise_plotcord, show.legend = FALSE, box.padding = 0.1) +
  theme_classic() + labs(x="", y="") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size = 8),
        plot.margin=unit(c(2, 2, 2, 2), "mm"))

plot(MEpaleturquoise_hubs)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (3) MEpaleturquoise g:Profiler enrichment ##############################################################################

library(gProfileR)

# important to reload the moduleColors so that we can get gene IDs
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')

query_background <- rownames(mmul_CON_vs_MIA_ALL)
#query_background <- query_background[which(query_background != '')]

MEpaleturquoise_ids <- names(moduleColors)[which(moduleColors == 'paleturquoise')]
MEpaleturquoise_kMEs <- kME$kMEpaleturquoise[which(rownames(kME) %in% MEpaleturquoise_ids)]
names(MEpaleturquoise_kMEs) <- rownames(kME)[which(rownames(kME) %in% MEpaleturquoise_ids)]
MEpaleturquoise_kMEs <- MEpaleturquoise_kMEs[order(MEpaleturquoise_kMEs, decreasing = TRUE)]

go = gprofiler(names(MEpaleturquoise_kMEs), organism="mmulatta", custom_bg = query_background,
               correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
               region_query = F,max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
               include_graph = F,src_filter = c("GO", "KEGG"))

go <- go[order(go$p.value, decreasing = FALSE),]
GO_summary1 <- data.frame('term' = go$term.name, 'p' = go$p.value, 'core_enrichment' = go$intersection, stringsAsFactors = FALSE)
GO_summary1 <- GO_summary1[c(1:5),]
GO_summary1$term <- c('neuron projection',
                      'regulation of plasma\nmembrane bounded cell\nprojection organization',
                      'gap junction hemi\nchannel activity',
                      'protein hexamerization',
                      'gap junction')
GO_summary1$term <- factor(GO_summary1$term, levels = GO_summary1$term)

MEpaleturquoise_GO <- ggplot(data = GO_summary1, aes(x = GO_summary1$term, y = -log10(GO_summary1$p))) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'paleturquoise') +
  scale_x_discrete(limits = rev(levels(droplevels(GO_summary1$term)))) +
  scale_y_continuous(breaks = c(0, 3, 6, 9), limits = c(0, 9)) +
  #facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 10, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(-2,5,5,5),"mm"))

plot(MEpaleturquoise_GO)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (4) Salmon module boxplot #######################################################################################

k <- which(colnames(MEs) == 'MEsalmon')

MEsalmon_boxplot <- ggplot(data.frame(MEs[,k]), aes(x=datMeta$MIA, y=MEs[,k])) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, size = 2, aes(color = datMeta$MIA)) +
  facet_grid(~ datMeta$Region) + 
  ggtitle(colnames(MEs)[k]) +
  theme_bw() +
  ylab(expression('Eigengene Expr')) +
  scale_color_discrete(labels = c("Control", "1st Trim MIA", "2nd Trim MIA")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(10,1,5,1),"mm"))


# plot(MEsalmon_boxplot)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (5) Salmon module hubs ##########################################################################################

source('./code/z_WGCNA_TOOLS.R')

datExpr <- mmul_txi.gene
datExpr <- datExpr[-which(duplicated(datExpr)),]
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr), method = 'TMM'), log = TRUE)
datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                           all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                           figureDir = './docs/trait_correlations/', figureName = 'mmul_gene_WGCNA')


library(igraph)

kME = signedKME(t(datExpr.norm.reg), MEs, corFnc = "bicor")

hubGenes <- rownames(kME)[order(kME[,which(colnames(kME) == 'kMEsalmon')], decreasing = T)[1:20]]
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


MEsalmon_hubs <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edges, size = 0.5, color = "grey") + 
  geom_point(aes(X1, X2), data = plotcord, color='salmon', size = 3) + 
  geom_text_repel(aes(x = X1, y = X2 + 0.1, label = gene, fontface = 'bold'), 
                  size = 2, data = plotcord, show.legend = FALSE, box.padding = 0.1) +
  theme_classic() + labs(x="", y="") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size = 8),
        plot.margin=unit(c(5, 5, 5, 5), "mm"))

plot(MEsalmon_hubs)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (6) MEsalmon g:Profiler enrichment ##############################################################################

library(gProfileR)

# important to reload the moduleColors so that we can get gene IDs
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')

query_background <- rownames(mmul_CON_vs_MIA_ALL)
#query_background <- query_background[which(query_background != '')]

MEsalmon_ids <- names(moduleColors)[which(moduleColors == 'salmon')]
MEsalmon_kMEs <- kME$kMEsalmon[which(rownames(kME) %in% MEsalmon_ids)]
names(MEsalmon_kMEs) <- rownames(kME)[which(rownames(kME) %in% MEsalmon_ids)]
MEsalmon_kMEs <- MEsalmon_kMEs[order(MEsalmon_kMEs, decreasing = TRUE)]

go = gprofiler(names(MEsalmon_kMEs), organism="mmulatta", custom_bg = query_background,
               correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
               region_query = F,max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
               include_graph = F,src_filter = c("GO", "KEGG"))

go <- go[order(go$p.value, decreasing = FALSE),]
GO_summary2 <- data.frame('term' = go$term.name, 'p' = go$p.value, 'core_enrichment' = go$intersection, stringsAsFactors = FALSE)
GO_summary2 <- GO_summary2[c(1:5),]
GO_summary2$term[5] <- 'oxidoreductase activity, acting on the CH-NH2 
                       \ngroup of donors, NAD or NADP as acceptor'
GO_summary2$term <- factor(GO_summary2$term, levels = GO_summary2$term)

MEsalmon_GO <- ggplot(data = GO_summary2, aes(x = GO_summary2$term, y = -log10(GO_summary2$p))) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'salmon') +
  scale_x_discrete(limits = rev(levels(droplevels(GO_summary2$term)))) +
  #facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  ylim(c(0, 6.5)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 10, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(-2,5,5,5),"mm"))

plot(MEsalmon_GO)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (7) Heatmaps of cell-type specificity for MEpaleturquoise ##########################################################

library(EWCE)

load('./working_data/CellTypeData_PEC_Lake_Adult_Human_nucSeq.rda')

library(reshape2)

head(ctd[[2]]$specificity)

MEpaleturquoise <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[match(MEpaleturquoise, mmul_annotEns87$ensembl_gene_id)]
MEpaleturquoise <- MEpaleturquoise[which(MEpaleturquoise != '')]

MEpaleturquoise_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% MEpaleturquoise),]
MEpaleturquoise_specificity <- MEpaleturquoise_specificity[order(apply(MEpaleturquoise_specificity, 1, max), decreasing = TRUE),]
# MEpaleturquoise_specificity <- MEpaleturquoise_specificity[1:50,]

MEpaleturquoise_specificity <- melt(MEpaleturquoise_specificity)
colnames(MEpaleturquoise_specificity) <- c('gene', 'cell', 'specificity')
MEpaleturquoise_specificity$specificity <- MEpaleturquoise_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

MEpaleturquoise_specificity_heatmap <- ggplot(MEpaleturquoise_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(MEpaleturquoise_specificity$specificity > 50, round(MEpaleturquoise_specificity$specificity, 1), '')), 
            size = 1.4, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  ggtitle('MEpaleturquoise') + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin=unit(c(2,0,5,0),"mm"))

plot(MEpaleturquoise_specificity_heatmap)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (8) turquoise module boxplot #######################################################################################

m <- which(colnames(MEs) == 'MEturquoise')

MEturquoise_boxplot <- ggplot(data.frame(HC_MEs[,m]), aes(x=HC_datMeta$MIA, y=HC_MEs[,m])) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 16, size = 2, aes(color = HC_datMeta$MIA)) +
  facet_grid(~ HC_datMeta$Region) + 
  ggtitle(colnames(HC_MEs)[m]) +
  theme_bw() +
  ylab(expression('Eigengene Expr')) +
  scale_color_discrete(labels = c("Control", "1st Trim MIA", "2nd Trim MIA")) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(2,1,2,1),"mm"))


# plot(MEturquoise_boxplot)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     MEturquoise_boxplot,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (9) turquoise module hubs ##########################################################################################

source('./code/z_WGCNA_TOOLS.R')

datExpr <- mmul_txi.gene
datExpr <- datExpr[-which(duplicated(datExpr)),]
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr), method = 'TMM'), log = TRUE)
datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                           all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                           figureDir = './docs/trait_correlations/', figureName = 'mmul_gene_WGCNA')


library(igraph)

kME = signedKME(t(datExpr.norm.reg), MEs, corFnc = "bicor")

hubGenes <- rownames(kME)[order(kME[,which(colnames(kME) == 'kMEturquoise')], decreasing = T)[1:20]]
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
MEturquoise_plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))


MEturquoise_hubs <- ggplot(data = MEturquoise_plotcord) + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data = edges, size = 0.5, color = "grey") + 
  geom_point(aes(X1, X2), data = MEturquoise_plotcord, color='turquoise', size = 3) + 
  geom_text_repel(aes(x = X1, y = X2 + 0.1, label = gene, fontface = 'bold'), 
                  size = 2, data = MEturquoise_plotcord, show.legend = FALSE, box.padding = 0.1) +
  theme_classic() + labs(x="", y="") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        text = element_text(size = 8),
        plot.margin=unit(c(2, 2, 2, 2), "mm"))

plot(MEturquoise_hubs)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     MEturquoise_boxplot, MEturquoise_hubs,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (10) MEturquoise g:Profiler enrichment ##############################################################################

library(gProfileR)

# important to reload the moduleColors so that we can get gene IDs
load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')

query_background <- rownames(mmul_CON_vs_MIA_ALL)
#query_background <- query_background[which(query_background != '')]

MEturquoise_ids <- names(moduleColors)[which(moduleColors == 'turquoise')]
MEturquoise_kMEs <- kME$kMEturquoise[which(rownames(kME) %in% MEturquoise_ids)]
names(MEturquoise_kMEs) <- rownames(kME)[which(rownames(kME) %in% MEturquoise_ids)]
MEturquoise_kMEs <- MEturquoise_kMEs[order(MEturquoise_kMEs, decreasing = TRUE)]

go = gprofiler(names(MEturquoise_kMEs), organism="mmulatta", custom_bg = query_background,
               correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
               region_query = F,max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
               include_graph = F,src_filter = c("GO", "KEGG"))

go <- go[order(go$p.value, decreasing = FALSE),]
GO_summary3 <- data.frame('term' = go$term.name, 'p' = go$p.value, 'core_enrichment' = go$intersection, stringsAsFactors = FALSE)
GO_summary3 <- GO_summary3[c(1:5),]
GO_summary3$term[4] <- 'structural constituent\nof myelin sheath'
GO_summary3$term <- factor(GO_summary3$term, levels = GO_summary3$term)

MEturquoise_GO <- ggplot(data = GO_summary3, aes(x = GO_summary3$term, y = -log10(GO_summary3$p))) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'turquoise') +
  scale_x_discrete(limits = rev(levels(droplevels(GO_summary3$term)))) +
  #facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  scale_y_continuous(breaks = c(0, 3, 6, 9), limits = c(0, 9)) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 10, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(-2,5,5,5),"mm"))

plot(MEturquoise_GO)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     MEturquoise_boxplot, MEturquoise_GO,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (11) MEturquoise cell-type enrichment ##############################################################################

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

MEturquoise = names(moduleColors[which(moduleColors == 'turquoise')])

background_adult = rownames(ctd[[level1]]$specificity)
length(background_adult)
background_adult = background_adult[which(background_adult %in% mmul_CON_vs_MIA_ALL$external_gene_name)]
length(background_adult)

MEturquoise = MEturquoise[which(MEturquoise %in% background_adult)]

MEturquoise_results = bootstrap.enrichment.test(sct_data=ctd, hits=MEturquoise, bg=background_adult,reps=reps,annotLevel=level1,genelistSpecies="human",sctSpecies="human")

ewce_results <- data.frame(MEturquoise_results$results[,c(1,3,5)], 'Module' = 'MEturquoise')

ewce_results$FDR <- p.adjust(ewce_results$p, method = 'fdr')
ewce_results$sd_from_mean <- ifelse(ewce_results$sd_from_mean > 0, ewce_results$sd_from_mean, 0)
print(ewce_results)

ewce_results$sig <- ifelse(ewce_results$FDR < 0.05, '*', '')

MEturquoise_cell_type_enrichment <- ggplot(data = ewce_results, aes(x = ewce_results$CellType, y = ewce_results$sd_from_mean, fill = ewce_results$Module)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_grid(cols = vars(Module)) +
  geom_text(aes(label=ewce_results$sig), vjust=0.5, color="red",
            position = position_dodge(0), size = 8) +
  scale_fill_manual(values = c('turquoise')) +
  ggtitle('Cell-type Enrichment\nfrom PsychEncode') +
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
        plot.margin=unit(c(5,2,5,2),"mm"))

plot(MEturquoise_cell_type_enrichment)

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     MEturquoise_boxplot, MEturquoise_GO, MEturquoise,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')



### (12) Combined WGCNA module DE plot #################################################################################

corfit <- duplicateCorrelation(t(MEs), mod, block = datMeta$Subject)
fit <- lmFit(t(MEs), mod, block = datMeta$Subject, correlation = corfit$consensus.correlation)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)
results <- eBayes(fit_contrasts)

# get the results of the linear model
tt_all <- topTable(results, number = Inf, coef = 1, sort.by = 'none')

tt_all <- tt_all[order(tt_all$logFC, decreasing = FALSE),]
tt_all$module <- factor(rownames(tt_all), levels = rownames(tt_all))

library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

combined_ME_barplot <- ggplot(tt_all, aes(x = module, y = logFC, fill = -log10(P.Value))) + 
  scale_fill_gradientn(colours = palette(100), limits = c(0, 1.5)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = c('*', rep('', 28))), vjust = 1.1, color="red",
            position = position_dodge(0), size=6) +
  ylim(c(-0.09, 0.09)) +
  ylab(expression('Effect Size (log'[2]*'(FC))')) +
  labs(fill = expression('-log'[10]*'(p)')) +
  ggtitle('Combined Module DE') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 0.5, hjust = 1, margin = margin(t = -15, r = 0, b = 20, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.margin = margin(t = -2, r = -1, b = -2, l = -2, unit = 'mm'))

save(MEpaleturquoise_boxplot, MEpaleturquoise_hubs, MEpaleturquoise_GO,
     MEsalmon_boxplot, MEsalmon_hubs, MEsalmon_GO, MEpaleturquoise_specificity_heatmap,
     MEturquoise_boxplot, MEturquoise_GO, MEturquoise, combined_ME_barplot,
     file = './working_data/plots/32_WGCNA_FIGURES.RData')


















