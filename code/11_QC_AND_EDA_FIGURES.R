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

########################################################################################################################
#
#   Pipeline:
#
#   (1) PCA plot colored by Region, MIA, and SeqBatch for Fig 1
#   (2) Trait correlations for Fig S1
#   (3) Pre-normalization batch effects for Fig S1
#   (4) Post-normalization batch effects for Fig S1
#   (5) Outlier Removal for Fig S1
#
########################################################################################################################

### (1) PCA plot colored by Region, MIA, and SeqBatch for Fig 1 ########################################################

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')
datMeta_outlier_removed <- datMeta[c(1:51),] #remove sample 40482_V1

mmul_txi.gene.norm <- cpm(calcNormFactors(DGEList(mmul_txi.gene), method = 'TMM'), log = TRUE)

ePCs <- prcomp(t(scale(t(mmul_txi.gene.norm), scale = TRUE)), center = FALSE)
ePCs <- data.frame(ePCs$rotation)

ePC_plot <- ggplot(ePCs, aes(x = PC1, y = PC2, color = datMeta_outlier_removed$Region, shape = datMeta_outlier_removed$SeqBatch, fill = datMeta_outlier_removed$MIA)) +
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

plot(ePC_plot)



### (2) Trait correlations for Fig S1 ##################################################################################

library(GGally)

pairsDat <- datMeta[,c(1:4,6:12)]
pairsDat <- cbind(pairsDat, ePCs[,c(1:5)])
colnames(pairsDat) <- c(colnames(pairsDat)[1:7],
                        substr(colnames(pairsDat)[8:11], 4, 9),
                        paste0('expr', colnames(pairsDat)[12:16]))

pairsDat$Subject <- as.numeric(pairsDat$Subject)
pairsDat$Region <- as.numeric(pairsDat$Region)
pairsDat$Group <- as.numeric(pairsDat$Group)
pairsDat$SeqBatch <- as.numeric(pairsDat$SeqBatch)
pairsDat$MIA <- as.numeric(pairsDat$MIA)
head(pairsDat)

trait_cor_plot <- ggpairs(pairsDat,
                          upper = list(continuous = wrap("cor", size = 1, color = 'black')),
                          lower = list(continuous = wrap("points", alpha = 1, size = 0.3))) +
  ggtitle("Trait Correlations -- |Pearson's rho|") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 3),
        axis.text.x = element_text(size = 4, angle = 270),
        axis.text.y = element_text(size = 4),
        plot.margin = unit(c(10, 10, 10, 10), 'mm'))

trait_cor_plot <- ggmatrix_gtable(trait_cor_plot)
plot(trait_cor_plot)

save(ePC_plot, trait_cor_plot, file = './working_data/plots/11_QC_AND_EDA_FIGURES.RData')



### (3) Pre-normalization batch effects for Fig S1 #####################################################################

library(reshape2)

melt_datExpr <- melt(log2(mmul_txi.gene + 1))
melt_datExpr$SeqBatch <- datMeta[match(melt_datExpr$Var2, rownames(datMeta)),]$SeqBatch

pre_norm_density <- ggplot(melt_datExpr, aes(x = value, color = Var2)) +
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

save(ePC_plot, trait_cor_plot, pre_norm_density, file = './working_data/plots/11_QC_AND_EDA_FIGURES.RData')



### (4) Post-normalization batch effects for Fig S1 ####################################################################

melt_datExpr <- melt(mmul_txi.gene.norm)
melt_datExpr$SeqBatch <- datMeta[match(melt_datExpr$Var2, rownames(datMeta)),]$SeqBatch

post_norm_density <- ggplot(melt_datExpr, aes(x = value, color = Var2)) +
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

# plot(post_norm_density)

save(ePC_plot, trait_cor_plot, pre_norm_density, post_norm_density, file = './working_data/plots/11_QC_AND_EDA_FIGURES.RData')



### (5) Outlier Removal for Fig S1 #####################################################################################

library(WGCNA)
library(ggrepel)

load('./raw_data/txi.kallisto.MMUL8v87.rda')
mmul_annotEns87 <- annotEns87
mmul_txi.gene <- txi.gene
mmul_txi.tx <- txi.tx

# import sample meta data and clean dataframe
datMeta <- read.csv(file = './raw_data/datMeta_20160601.csv', header = TRUE, sep = ',')
rownames(datMeta) <- datMeta$Sample
datMeta <- datMeta[,c(5:11,16:21)]
colnames(datMeta) <- c(colnames(datMeta)[1:8], 'MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4', 'MG_SeqDepth')
datMeta$Subject <- as.factor(datMeta$Subject)
datMeta$Region <- factor(datMeta$Region, levels = c('DLPFC', 'ACC', 'HC', 'V1'))
datMeta$Group <- factor(datMeta$Group, levels = c('con', 'poly1', 'poly2'))
datMeta$SeqBatch <- factor(datMeta$SeqBatch, levels = c('1', '2'))
datMeta$MIA <- factor(datMeta$MIA, levels = c('CON', 'MIA'))

genes_to_keep <- filterByExpr(DGEList(mmul_txi.gene$counts)) #need to mention this in the methods section
table(genes_to_keep)
mmul_txi.gene <- mmul_txi.gene$counts[genes_to_keep,]

mmul_txi.gene.norm <- cpm(calcNormFactors(DGEList(mmul_txi.gene), method = 'TMM'), log = TRUE)


# get connectivity z-scores for the normalized data
normadj <- (0.5 + 0.5*bicor(mmul_txi.gene.norm)^2)
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

# plot the connectivity z-score for the normalized data and set a cut-off at 3 STD
plot(1:length(Z.C),Z.C,main=paste0("Outlier Plot of"),xlab = "\nSamples\n",ylab="\nConnectivity Z Score\n", ylim = c(-4, 1), col = datMeta$SeqBatch)
legend("topright", title = 'SeqBatch', legend = levels(datMeta$SeqBatch), fill = rev(as.numeric(as.factor(levels(datMeta$SeqBatch)))))
abline(h = -3, col="red")

connectivity <- data.frame(sample_num = 1:length(Z.C), z_score = Z.C, batch = ifelse(datMeta$SeqBatch == 1, 'Batch 1', 'Batch 2'), name = rownames(datMeta))

outlier_plot <- ggplot(connectivity, aes(x = sample_num, y = z_score, color = batch)) +
  geom_point() +
  geom_hline(yintercept = -3, color = 'blue') + 
  scale_color_manual(values = c('black', 'red')) +
  ggtitle('Outlier Removal') +
  xlab('Sample #') +
  ylab('Connectivity Z-Score') +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(color = c('black', 'red'), shape = c(16, 16)))) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(size = 8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, 'mm'),
        legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.direction = 'vertical',
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  geom_text_repel(data = data.frame(connectivity[52,]), aes(x = sample_num, y = z_score, label = name), size = 2)

# plot(outlier_plot)

save(ePC_plot, trait_cor_plot, pre_norm_density, post_norm_density, outlier_plot, file = './working_data/plots/11_QC_AND_EDA_FIGURES.RData')















