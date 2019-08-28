########################################################################################################################
#
#   41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.R 
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
load('./working_data/DGE_result_topTables.RData')
load('./working_data/TE_result_topTables.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Plot TE PCA
#   (2) Plot TE expression density pre-normalization
#   (3) Plot TE expression density post-normalization
#   (4) Plot global differentially expressed TEs
#   (5) Plot PGBD3 expression in ACC
#   (6) Plot PIWIL2 expression across all regions
#
########################################################################################################################

### (1) Plot TE PCA ####################################################################################################

rep_datExpr <- read.csv(file = './raw_data/TransposableElements/RepEnrich2/TE_fraction_counts.csv', header = TRUE, sep = '\t')
rownames(rep_datExpr) <- rep_datExpr$TE
rep_datExpr <- rep_datExpr[,c(2:53)]
colnames(rep_datExpr) <- substring(colnames(rep_datExpr), 2)

# remove outlier V1 sample from all dataframes
rep_datExpr <- rep_datExpr[,which(colnames(rep_datExpr) != '40482_V1')]

# remove low expresssed TEs
genes_to_keep <- filterByExpr(DGEList(rep_datExpr)) #need to mention this in the methods section
table(genes_to_keep)
rep_datExpr <- rep_datExpr[genes_to_keep,] # removed 72 low expressed TEs

rep_datExpr.norm <- cpm(calcNormFactors(DGEList(rep_datExpr), method = 'TMM'), log = TRUE)

# plot first two expression PCs
TE_ePCs <- prcomp(t(scale(t(rep_datExpr.norm), scale = TRUE)), center = FALSE)
TE_ePCs <- data.frame(TE_ePCs$rotation)
TE_ePC_plot <- ggplot(TE_ePCs, aes(x = PC1, y = PC2, color = datMeta$Region, shape = datMeta$SeqBatch, fill = datMeta$MIA)) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(22, 21), labels = c('Batch 1', 'Batch 2')) +
  scale_fill_manual(values = c('white', 'grey'), labels = c('Batch 1', 'Batch 2')) +
  ggtitle('TE Quality Control') +
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

#plot(TE_ePC_plot)

save(TE_ePC_plot,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (2) Plot TE expression density pre-normalization ###################################################################

library(reshape2)

melt_datExpr <- melt(log2(rep_datExpr + 1))
melt_datExpr$SeqBatch <- datMeta[match(melt_datExpr$variable, rownames(datMeta)),]$SeqBatch

TE_pre_norm_density <- ggplot(melt_datExpr, aes(x = value, color = variable)) +
  geom_line(stat = 'density', size = 0.02) +
  geom_point(data = data.frame(point = c(10, 10), seqbatch = c('Batch 1', 'Batch 2')), aes(x = point[1], y = point[2], shape = seqbatch, color = 'white'), inherit.aes = FALSE) +
  scale_color_manual(values = c(datMeta$SeqBatch, 'white'), labels = c(datMeta$SeqBatch, 'white'), guide = 'none') +
  xlim(c(-1, 20)) +
  ylim(c(0, 0.2)) + 
  theme_bw() +
  ggtitle('Pre-Normalization') +
  xlab(expression('log'[2]*'(counts+1)')) +
  ylab('Density') +
  # guides(shape = guide_legend(override.aes = list(shape = c(16, 16), color = c('black', 'red')))) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_blank(),
        # legend.text = element_text(size = 6),
        legend.position = 'none',
        legend.text = element_blank(),legend.key.size = unit(0.5, 'mm'),
        legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.direction = 'vertical',
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

# plot(pre_norm_density)

save(TE_ePC_plot, TE_pre_norm_density,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (3) Plot TE expression density post-normalization ##################################################################

melt_datExpr <- melt(rep_datExpr.norm)
melt_datExpr$SeqBatch <- datMeta[match(melt_datExpr$Var2, rownames(datMeta)),]$SeqBatch

TE_post_norm_density <- ggplot(melt_datExpr, aes(x = value, color = Var2)) +
  geom_line(stat = 'density', size = 0.02) +
  geom_point(data = data.frame(point = c(10, 10), seqbatch = c('Batch 1', 'Batch 2')), aes(x = point[1], y = point[2], shape = seqbatch, color = 'white'), inherit.aes = FALSE) +
  scale_color_manual(values = c(datMeta$SeqBatch, 'white'), labels = c(datMeta$SeqBatch, 'white'), guide = 'none') +
  xlim(c(-5, 20)) +
  ylim(c(0, 0.2)) + 
  theme_bw() +
  ggtitle('Post-Normalization') +
  xlab(expression('log'[2]*'(counts+1)')) +
  ylab('Density') +
  # guides(shape = guide_legend(override.aes = list(shape = c(16, 16), color = c('black', 'red')))) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_blank(),
        # legend.text = element_text(size = 6),
        legend.position = 'none',
        legend.text = element_blank(),legend.key.size = unit(0.5, 'mm'),
        legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.direction = 'vertical',
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

# plot(TE_post_norm_density)

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (4) Plot global differentially expressed TEs #######################################################################

mmul_CON_vs_MIA_ALL <- TE_result_topTables$RepEnrich_CON_vs_MIA_ALL

# use and ifelse statement to assign the correct color to genes in the volcano plot based on their P-Values
to_color_all <- c(ifelse(mmul_CON_vs_MIA_ALL$FDR < 0.1, 1, ifelse(mmul_CON_vs_MIA_ALL$P.Value < 0.005 & mmul_CON_vs_MIA_ALL$FDR >= 0.1 ,2, 3)))
to_color_table_all <- table(to_color_all)

# copy ensembl ID to external gene name if it is missing
mmul_CON_vs_MIA_ALL$external_gene_name <- ifelse(mmul_CON_vs_MIA_ALL$external_gene_name == '', '', mmul_CON_vs_MIA_ALL$external_gene_name)

mmul_all <- mmul_CON_vs_MIA_ALL

# generate Limma-voom volcano plot for each region and time point
TE_combined_volcano_plot <- ggplot(as.data.frame(mmul_all), aes(x=mmul_all$logFC, y=-log10(mmul_all$P.Value), color = as.factor(to_color_all), label=mmul_all$external_gene_name)) + 
  geom_point(alpha = 1, size = 2) +  scale_color_manual(labels = c('1' = paste0('FDR < 0.1 (', to_color_table_all[names(to_color_table_all) == 1], ' TEs)'), 
                                                                   '2' = paste0('P-val < 0.005 (', to_color_table_all[names(to_color_table_all) == 2], ' TEs)'), 
                                                                   '3' = paste0('n.s. (', to_color_table_all[names(to_color_table_all) == 3], ' TEs)')),
                                                        values=c('1' = "red", '2' = "orange", '3' = "grey60"))  + 
  guides(fill=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) + 
  ggtitle(paste0('Global Transposable Element DE')) + xlab(expression('Effect Size (log'[2]*'(FC))')) + ylab(expression('\nSignificance (-log'[10]*'(P-Value))')) + labs(color = '') +
  geom_text_repel(data = mmul_all[1:8,], aes(x=mmul_all[1:8,]$logFC, y=-log10(mmul_all[1:8,]$P.Value), label = mmul_all[1:8,]$external_gene_name), inherit.aes = FALSE,
                  size = 3, color="black", alpha=1, box.padding = unit(0.45, "lines")) +
  theme_classic() + theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                          axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                          axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                          legend.text=element_text(size = 6),
                          axis.text.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),
                          #legend.justification = c(1, 1),
                          #legend.position = c(0.2, 1.08),
                          legend.direction = 'horizontal',
                          legend.position = 'bottom',
                          plot.margin = unit(c(0, 0, 0, 0), 'mm'),
                          legend.margin = margin(t = -2, r = 0, b = 0, l = 0, unit = 'mm')) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 1)

# plot(TE_combined_volcano_plot)

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     TE_combined_volcano_plot,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (5) Plot PGBD3 expression in ACC ###################################################################################


library(scales)

tt.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
DLPFC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_V1


datExpr <- mmul_txi.gene
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr)), method = 'TMM', log = TRUE)
to_keep <- which(datMeta$Region == 'ACC')
ACC_datMeta <- datMeta[to_keep,]
ACC_datExpr.norm <- datExpr.norm[,to_keep]

TE_genes <- c('PGBD3')
TE_gene_ids <- mmul_annotEns87$ensembl_gene_id[match(TE_genes, mmul_annotEns87$hsapiens_homolog_associated_gene_name)]


PGBD3 <- ggplot(data.frame(ACC_datExpr.norm[which(rownames(ACC_datExpr.norm) == TE_gene_ids[1]),]), aes(x=ACC_datMeta$MIA, y=ACC_datExpr.norm[which(rownames(ACC_datExpr.norm) == TE_gene_ids[1]),])) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75), aes(color = ACC_datMeta$MIA)) +
  facet_grid(. ~ ACC_datMeta$Region) + 
  ggtitle(TE_genes[1], subtitle = paste0('ACC: p=', round(ACC.fullModel$P.Value[which(ACC.fullModel$external_gene_name == TE_genes[1])], 4))) +
  theme_bw() +
  ylab(expression('Expr (log '[2]*' (CPM))')) +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 6),
        legend.position = 'none',
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     TE_combined_volcano_plot, PGBD3,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (6) Plot PIWIL2 expression across all regions ######################################################################

library(scales)

tt.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
DLPFC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_V1

datExpr <- mmul_txi.gene
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr)), method = 'TMM', log = TRUE)

genes <- c('PIWIL2')
gene_ids <- mmul_annotEns87$ensembl_gene_id[match(genes, mmul_annotEns87$hsapiens_homolog_associated_gene_name)]


PIWIL2 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[1]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[1]),])) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6, color = 'black') +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75), aes(color = datMeta$MIA)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[1]) +
  theme_bw() +
  ylab(expression('Expr (log '[2]*' (CPM))')) +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 6),
        legend.position = 'none')

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     TE_combined_volcano_plot, PGBD3, PIWIL2,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (7) Regional TE DE volcano plot ####################################################################################

mmul_CON_vs_MIA_DLPFC <- TE_result_topTables$RepEnrich_CON_vs_MIA_DLPFC
mmul_CON_vs_MIA_DLPFC['Region'] <- 1

mmul_CON_vs_MIA_ACC <- TE_result_topTables$RepEnrich_CON_vs_MIA_ACC
mmul_CON_vs_MIA_ACC['Region'] <- 2

mmul_CON_vs_MIA_HC <- TE_result_topTables$RepEnrich_CON_vs_MIA_HC
mmul_CON_vs_MIA_HC['Region'] <- 3

mmul_CON_vs_MIA_V1 <- TE_result_topTables$RepEnrich_CON_vs_MIA_V1
mmul_CON_vs_MIA_V1['Region'] <- 4

mmul_CON_vs_MIA_Combined <- rbind(mmul_CON_vs_MIA_DLPFC,
                                  mmul_CON_vs_MIA_ACC,
                                  mmul_CON_vs_MIA_HC,
                                  mmul_CON_vs_MIA_V1)


mmul_CON_vs_MIA_Combined['all_FDR'] <- fdrtool::fdrtool(mmul_CON_vs_MIA_Combined$t, plot = TRUE)$qval

top_mmul_CON_vs_MIA_Combined <- rbind(mmul_CON_vs_MIA_DLPFC[1:9,],
                                      mmul_CON_vs_MIA_ACC[1:6,],
                                      mmul_CON_vs_MIA_HC[1:12,])

# use and ifelse statement to assign the correct color to genes in the volcano plot based on their P-Values
to_color <- c(ifelse(mmul_CON_vs_MIA_Combined$all_FDR < 0.1, 1, ifelse(mmul_CON_vs_MIA_Combined$P.Value < 0.005 & mmul_CON_vs_MIA_Combined$all_FDR >= 0.1 ,2, 3)))
to_color_table <- table(to_color)

region.labels <- c('DLPFC', 'ACC', 'HC', 'V1')
names(region.labels) <- c(1, 2, 3, 4)

TE_regional_volcano_plot <- ggplot(as.data.frame(mmul_CON_vs_MIA_Combined), aes(x=mmul_CON_vs_MIA_Combined$logFC, y=-log10(mmul_CON_vs_MIA_Combined$P.Value), color = as.factor(to_color), label=mmul_CON_vs_MIA_Combined$external_gene_name)) + 
  facet_wrap(~ Region, nrow = 1, labeller = labeller(Region = region.labels)) +
  geom_point(alpha=1, size = 2) +  scale_color_manual(labels = c('1' = paste0('FDR < 0.1'), 
                                                                 '2' = paste0('P-val < 0.005'), 
                                                                 '3' = paste0('n.s.')),
                                                      values=c('1' = "red", '2' = "orange", '3' = "grey60"))  + 
  guides(fill=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) + 
  ggtitle(paste0('Regional Transposable Element DE')) + xlab(expression('Effect Size (log'[2]*'(FC))')) + ylab(expression('\nSignificance (-log'[10]*'(P-Value))')) + labs(color = '') +
  geom_text_repel(data = top_mmul_CON_vs_MIA_Combined, aes(x=top_mmul_CON_vs_MIA_Combined$logFC, y=-log10(top_mmul_CON_vs_MIA_Combined$P.Value), label = top_mmul_CON_vs_MIA_Combined$external_gene_name), inherit.aes = FALSE,
                  size = 3, color="black", alpha=1, box.padding = unit(0.65, "lines")) +
  theme_bw() + theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
                     axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                     legend.text = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     strip.text.x = element_text(size = 10),
                     legend.margin = margin(0, 0, 0, 0),
                     legend.box.margin = margin(-10, -10, -10, -10),
                     legend.position = 'bottom',
                     legend.direction = 'horizontal')

# plot(TE_regional_volcano_plot)

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     TE_combined_volcano_plot, PGBD3, PIWIL2, TE_combined_volcano_plot,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')



### (5) Plot PGBD3 expression in ACC ###################################################################################


library(scales)

tt.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
DLPFC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_V1


datExpr <- mmul_txi.gene
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr)), method = 'TMM', log = TRUE)
to_keep <- which(datMeta$Region == 'HC')
HC_datMeta <- datMeta[to_keep,]
HC_datExpr.norm <- datExpr.norm[,to_keep]

HC_TE_genes <- c('PGBD5')
HC_TE_gene_ids <- mmul_annotEns87$ensembl_gene_id[match(HC_TE_genes, mmul_annotEns87$hsapiens_homolog_associated_gene_name)]


PGBD5 <- ggplot(data.frame(HC_datExpr.norm[which(rownames(HC_datExpr.norm) == HC_TE_gene_ids[1]),]), aes(x=HC_datMeta$MIA, y=HC_datExpr.norm[which(rownames(HC_datExpr.norm) == HC_TE_gene_ids[1]),])) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75), aes(color = HC_datMeta$MIA)) +
  facet_grid(. ~ HC_datMeta$Region) + 
  ggtitle(HC_TE_genes[1], subtitle = paste0('HC: p=', round(HC.fullModel$P.Value[which(HC.fullModel$external_gene_name == HC_TE_genes[1])], 4))) +
  theme_bw() +
  ylab(expression('Expr (log '[2]*' (CPM))')) +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 6),
        legend.position = 'none',
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

save(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density,
     TE_combined_volcano_plot, PGBD3, PGBD5,
     file = './working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')


