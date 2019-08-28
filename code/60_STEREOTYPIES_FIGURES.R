########################################################################################################################
#
#   60_STEREOTYPIES_FIGURES.R 
#
#   Determine cross disorder overlap between differential gene expression in primate MIA model
#   and neurodevelopmental disorders autism and schizophrenia
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

stereo <- read.csv('./raw_data/stereoAMA27.csv', header = TRUE, sep = ',')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Perform t-test on CON vs MIA total stereotypies
#   (2) Plot correlations between stereotypies and module eigengenes
#   (3) Plot TE expression density post-normalization
#   (4) RRHO on DLPFC DE genes vs disease
#
########################################################################################################################

### (1) Perform t-test on CON vs MIA total stereotypies ################################################################

stereo <- stereo[which(stereo$Sex == 'M'),]

result <- t.test(log2(stereo$TotalStereo_Sum[which(stereo$Poly_vs_Con == 'C')]),
                 log2(stereo$TotalStereo_Sum[which(stereo$Poly_vs_Con == 'P')]),
                 alternative = 'two.sided')

p_val <- result$p.value

stereo_plot <- ggplot(stereo, aes(x=stereo$Poly_vs_Con, y = log2(stereo$TotalStereo_Sum))) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 2) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75), aes(color = stereo$Poly_vs_Con)) +
  ggtitle('Total Stereotypies', subtitle = paste0('p: ', round(p_val, 3))) +
  scale_x_discrete(labels = c('Ctrl', 'MIA')) +
  theme_bw() +
  ylab('log(Total Stereotypies)') +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 8,  margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 6),
        axis.title.x=element_blank(),
        legend.position = 'none', 
        plot.margin = unit(c(3, 3, 3, 3), 'mm'))

plot(stereo_plot)

save(stereo_plot,
     file = './working_data/plots/60_STEREOTYPIES_FIGURES.RData')



### (2) Plot correlations between stereotypies and module eigengenes ###################################################

load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')

DLPFC_MEs <- MEs[which(datMeta$Region == 'DLPFC'),]
DLPFC_cor_mat <- cor(DLPFC_MEs, log2(stereo[,10]))

ACC_MEs <- MEs[which(datMeta$Region == 'ACC'),]
ACC_cor_mat <- cor(ACC_MEs, log2(stereo[,10]))

HC_MEs <- MEs[which(datMeta$Region == 'HC'),]
HC_cor_mat <- cor(HC_MEs, log2(stereo[,10]))

V1_MEs <- MEs[which(datMeta$Region == 'V1'),]
V1_cor_mat <- cor(V1_MEs, log2(stereo[c(1:12),10]))

colnames(DLPFC_cor_mat) <- 'DLPFC'
colnames(ACC_cor_mat) <- 'ACC'
colnames(HC_cor_mat) <- 'HC'
colnames(V1_cor_mat) <- 'V1'

cor_mat <- cbind(DLPFC_cor_mat, ACC_cor_mat, HC_cor_mat, V1_cor_mat)
melt_cor_mat <- melt(cor_mat)

library(RColorBrewer)
palette <- colorRampPalette(rev(brewer.pal(11, 'RdBu')), space='Lab')

stereo_cor_mat <- ggplot(melt_cor_mat, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(melt_cor_mat$value > 0.5, round(melt_cor_mat$value, 2), '')), 
            size = 2, color = 'black', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(-1, 1)) +
  ggtitle('Module Eigengene Stereotypies Correlation') +
  labs(fill = 'Correlation (R)') +
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
        legend.title = element_text(size = 8),
        plot.margin=unit(c(2,5,5,2),"mm"),
        legend.box.margin=margin(-10,5,-10,-10))

plot(stereo_cor_mat)

save(stereo_plot, stereo_cor_mat,
     file = './working_data/plots/60_STEREOTYPIES_FIGURES.RData')



### (3) Plot stereotypies vs MEdarkorange in V1 ########################################################################

MEdarkorange <- cbind(V1_MEs$MEdarkorange, stereo[c(1:12),c(2,10)])
MEdarkorange$total_stereo <- log2(MEdarkorange$total_stereo)
rownames(MEdarkorange) <- rownames(V1_MEs)
colnames(MEdarkorange) <- c('ME', 'Group', 'total_stereo')
MEdarkorange <- data.frame(MEdarkorange)

MEdarkorange_stereo_plot <- ggplot(MEdarkorange, aes(x = ME, y = total_stereo)) +
  geom_point(size = 2, aes(color = Group)) + 
  geom_smooth(method='lm', se = FALSE) +
  scale_color_manual(values = c('black', 'red'), labels = c('Ctrl', 'MIA')) +
  ggtitle('MEdarkorange vs. Total Stereotypies') +
  labs(color = 'Group') + 
  xlab('Eigengene Expr') +
  ylab('Total # Stereotypies') +
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

plot(MEdarkorange_stereo_plot)

save(stereo_plot, stereo_cor_mat, MEdarkorange_stereo_plot,
     file = './working_data/plots/60_STEREOTYPIES_FIGURES.RData')


