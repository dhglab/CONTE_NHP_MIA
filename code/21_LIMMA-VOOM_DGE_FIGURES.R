#####################################################################################################################################
#
#   21_CONTE_NHP_MIA_DGE_FIGURES.R 
#
#   Generate figures for differential gene expression on 4 y/o non-human 
#   primate RNA-seq samples from DLPFC, ACC, HC, and V1 (Part 1)
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
#####################################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(gridExtra)
library(ggplot2)
library(ggrepel)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')
load('./working_data/DGE_result_topTables.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Volcano plot of global differential gene expression
#   (2) g:Profiler enrichment for global up and down genes with p < 0.005
#   (3) Cell-type enrichment from PsychEncode for global up and down genes with p < 0.005
#   (4) Volcano plots of poly1 vs poly2 differential gene expression
#   (5) Global logFC correlation of top poly1 vs poly2 DE genes
#   (6) g:Profiler enrichment for top poly1 vs poly2 specific DE genes
#
########################################################################################################################

### (1) Volcano plot of global differential gene expression ############################################################

mmul_CON_vs_MIA_ALL <- DGE_result_topTables$mmul_CON_vs_MIA_ALL

# use and ifelse statement to assign the correct color to genes in the volcano plot based on their P-Values
to_color_all <- c(ifelse(mmul_CON_vs_MIA_ALL$FDR < 0.1, 1, ifelse(mmul_CON_vs_MIA_ALL$P.Value < 0.005 & mmul_CON_vs_MIA_ALL$FDR >= 0.1 ,2, 3)))
to_color_table_all <- table(to_color_all)

# copy ensembl ID to external gene name if it is missing
mmul_CON_vs_MIA_ALL$external_gene_name <- ifelse(mmul_CON_vs_MIA_ALL$external_gene_name == '', '', mmul_CON_vs_MIA_ALL$external_gene_name)

mmul_all <- mmul_CON_vs_MIA_ALL

# generate Limma-voom volcano plot for each region and time point
combined_volcano_plot <- ggplot(as.data.frame(mmul_all), aes(x=mmul_all$logFC, y=-log10(mmul_all$P.Value), color = as.factor(to_color_all), label=mmul_all$external_gene_name)) + 
  geom_point(alpha = 1, size = 2) +  scale_color_manual(labels = c('1' = paste0('FDR < 0.1 (', to_color_table_all[names(to_color_table_all) == 1], ' genes)'), 
                                                         '2' = paste0('P-val < 0.005 (', to_color_table_all[names(to_color_table_all) == 2], ' genes)'), 
                                                         '3' = paste0('n.s. (', to_color_table_all[names(to_color_table_all) == 3], ' genes)')),
                                              values=c('1' = "red", '2' = "orange", '3' = "grey60"))  + 
  guides(fill=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) + 
  ggtitle(paste0('Global Differential Gene Expression')) + xlab(expression('Effect Size (log'[2]*'(FC))')) + ylab(expression('\nSignificance (-log'[10]*'(P-Value))')) + labs(color = '') +
  geom_text_repel(data = mmul_all[1:50,], aes(x=mmul_all[1:50,]$logFC, y=-log10(mmul_all[1:50,]$P.Value), label = mmul_all[1:50,]$external_gene_name), inherit.aes = FALSE,
                 size = 3, color="black", alpha=1, box.padding = unit(0.45, "lines")) +
  theme_classic() + theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
                     axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                     legend.text=element_text(size = 8),
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     legend.justification = c(1, 1),
                     legend.position = c(1, 1.08),
                     legend.direction = 'vertical') +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 1)
  
plot(combined_volcano_plot)

save(combined_volcano_plot, file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')



### (2) g:Profiler enrichment for global up and down genes with p < 0.005 ##############################################

library(gProfileR)

# order dataset by -log10(signed P-Value) to use for g:Profiler enrichment
mmul_CON_vs_MIA_ALL$signed_p <- -log10(mmul_CON_vs_MIA_ALL$P.Value) * sign(mmul_CON_vs_MIA_ALL$logFC)
mmul_CON_vs_MIA_ALL_ordered <- mmul_CON_vs_MIA_ALL[order(mmul_CON_vs_MIA_ALL$signed_p, decreasing = TRUE),]

# get the top Up and Down regulated genes from the combined differential gene expression
combined_UP = rownames(mmul_CON_vs_MIA_ALL_ordered)[which(mmul_CON_vs_MIA_ALL_ordered$P.Value < 0.05 & mmul_CON_vs_MIA_ALL_ordered$logFC > 0)]
combined_DOWN = rev(rownames(mmul_CON_vs_MIA_ALL_ordered)[which(mmul_CON_vs_MIA_ALL_ordered$P.Value < 0.05 & mmul_CON_vs_MIA_ALL_ordered$logFC < 0)]) # need to take the reverse of this
                                                                                                                                                      # so genes with lowest signed-Pval 
                                                                                                                                                      # are first!!!
# Remove genes from gene list that lack a human homolog
combined_UP <- combined_UP[which(combined_UP != '')]
combined_DOWN <- combined_DOWN[which(combined_DOWN != '')]

query_background <- rownames(mmul_CON_vs_MIA_ALL_ordered)

query_list <- list(combined_UP, combined_DOWN)
query_region_list <- list('Combined', 'Combined')
query_direction_list <- list('UP', 'DOWN')

GO_combined <- list()
top_GO_combined <- list()

mapply(function(query, query_region, query_direction) {
  
  # Perform g:Profiler enrichment based off of an ordered query for mmulatta, max set size = 1000
  go = gprofiler(query, organism="mmulatta", custom_bg = query_background,
                 correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
                 region_query = F, max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
                 include_graph = F,src_filter = c("GO", "KEGG"))
  GO_combined <<- rbind(GO_combined, data.frame('term' = go$term.name, 'p' = go$p.value, 'Region' = query_region, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  go = go[order(go$p.value)[1:min(5,nrow(go))],]
  top_GO_combined <<- rbind(top_GO_combined, data.frame('term' = go$term.name, 'p' = go$p.value, 'Region' = query_region, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  
  
}, query_list, query_region_list, query_direction_list)

GO_combined <- data.frame(GO_combined)
head(GO_combined)
hist(GO_combined$p)
GO_combined <- GO_combined[order(GO_combined$p, decreasing = FALSE),]
head(GO_combined)

top_GO_combined <- data.frame(top_GO_combined)
top_GO_combined

GO_combined <- top_GO_combined

library(ggplot2)
library(lemon)

GO_combined_UP <- GO_combined[which(GO_combined$Direction == 'UP'),]
GO_combined_UP$term <- c('membrane fusion involved\nin viral entry into host cell',
                        'ribonuclease A activity',
                        'bile acid binding',
                        'nucleosome',
                        'cell-cell junction')

p1 <- ggplot(data = GO_combined_UP, aes(x = GO_combined_UP$term, y = -log10(GO_combined_UP$p), fill = GO_combined_UP$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_combined_UP$term)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_blank(),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10))

# plot(p1)


GO_combined_DOWN <- GO_combined[which(GO_combined$Direction == 'DOWN'),]
GO_combined_DOWN$term <- c('meiosis I cell cycle process',
                          'snoRNA binding',
                          '90S preribosome',
                          'germ-line stem cell\npopulation maintenance',
                          'perinucleolar chromocenter')

p2 <- ggplot(data = GO_combined_DOWN, aes(x = GO_combined_DOWN$term, y = -log10(GO_combined_DOWN$p), fill = GO_combined_DOWN$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_combined_DOWN$term)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0, 3)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ylab(expression('\n-log'[10]*'(FDR)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(-2,5,5,5),"mm"))

# plot(p2)


library(cowplot)

combined_GO_enrichment_plot <- plot_grid(p1, p2, ncol = 1, align = 'v')

save(combined_volcano_plot, combined_GO_enrichment_plot, file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')



### (3) Cell-type enrichment from Lake et al for global up and down genes with p < 0.005 ###############################

library(EWCE)

load('./working_data/CellTypeData_PEC_Lake_Adult_Human_nucSeq.rda')

set.seed(12211991)
reps = 10000

# get the top Up and Down regulated genes from the combined differential gene expression
combined_UP = mmul_CON_vs_MIA_ALL_ordered$external_gene_name[which(mmul_CON_vs_MIA_ALL_ordered$P.Value < 0.05 & mmul_CON_vs_MIA_ALL_ordered$logFC > 0)]
combined_DOWN = mmul_CON_vs_MIA_ALL_ordered$external_gene_name[which(mmul_CON_vs_MIA_ALL_ordered$P.Value < 0.05 & mmul_CON_vs_MIA_ALL_ordered$logFC < 0)] # need to take the reverse of this
                                                                                                                                                          # so genes with lowest signed-Pval 
                                                                                                                                                          # are first!!!
# Remove genes from gene list that lack a human homolog
combined_UP <- combined_UP[which(combined_UP != '')]
combined_DOWN <- combined_DOWN[which(combined_DOWN != '')]


background_adult = rownames(ctd[[2]]$specificity)
length(background_adult)
background_adult = background_adult[which(background_adult %in% mmul_CON_vs_MIA_ALL$external_gene_name)]
length(background_adult) # final set has 12228 genes

# remove genes that are not in the background for accurate results
combined_UP_EWCE = combined_UP[which(combined_UP %in% background_adult)]
combined_DOWN_EWCE = combined_DOWN[which(combined_DOWN %in% background_adult)]

# perform EWCE enrichment
combined_UP_EWCE_results = bootstrap.enrichment.test(sct_data = ctd, hits = combined_UP_EWCE, bg = background_adult, reps = reps, annotLevel = 2, genelistSpecies = 'human', sctSpecies = 'human')
combined_DOWN_EWCE_results = bootstrap.enrichment.test(sct_data = ctd, hits = combined_DOWN_EWCE, bg = background_adult, reps = reps, annotLevel = 2, genelistSpecies = 'human', sctSpecies = 'human')

ewce_results <- rbind(data.frame(combined_UP_EWCE_results$results[,c(1,3,5)], 'Region' = '1', 'Direction' = '1'),
                      data.frame(combined_DOWN_EWCE_results$results[,c(1,3,5)], 'Region' = '1', 'Direction' = '2'))

ewce_results$FDR <- p.adjust(ewce_results$p, method = 'fdr')
ewce_results$sd_from_mean <- ifelse(ewce_results$sd_from_mean > 0, ewce_results$sd_from_mean, 0)
print(ewce_results)

ewce_results$Direction <- factor(ewce_results$Direction, labels = c('UP', 'DOWN'))
ewce_results$Region <- factor(ewce_results$Region, labels = c('Combined'))
ewce_results$sig <- ifelse(ewce_results$FDR < 0.05, '*', '')

combined_EWCE_enrichment_plot <- ggplot(data = ewce_results, aes(x = ewce_results$CellType, y = ewce_results$sd_from_mean, fill = ewce_results$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  facet_wrap(~Direction, nrow = 1, strip.position = 'top') +
  scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 10)) +
  geom_text(aes(label = ewce_results$sig), vjust=0.5, color="red",
            position = position_dodge(0), size=8) +
  ggtitle('Cell-type Enrichment\nfrom PsychEncode') +
  ylab('Stdev from Mean') +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(10,5,10,5),"mm"))

# plot(combined_EWCE_enrichment_plot)

save(combined_volcano_plot, combined_GO_enrichment_plot, combined_EWCE_enrichment_plot, 
     file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')



### (4) Volcano plots of poly1 vs poly2 differential gene expression ###################################################

poly1_ALL <- DGE_result_topTables$mmul_con_vs_poly1
poly1_ALL['analysis'] <- '1st Trim MIA'
poly2_ALL <- DGE_result_topTables$mmul_con_vs_poly2
poly2_ALL['analysis'] <- '2nd Trim MIA'

poly_combined <- rbind(poly1_ALL,
                       poly2_ALL)

poly_combined['all_FDR'] <- fdrtool::fdrtool(poly_combined$t, plot = TRUE)$qval

top_poly_combined <- rbind(poly1_ALL[1:30,],
                           poly2_ALL[1:30,])

# use and ifelse statement to assign the correct color to genes in the volcano plot based on their P-Values
to_color <- c(ifelse(poly_combined$all_FDR < 0.1, 1, ifelse(poly_combined$P.Value < 0.005 & poly_combined$all_FDR >= 0.1 ,2, 3)))
to_color_table <- table(to_color)


poly_volcano_plot <- ggplot(as.data.frame(poly_combined), aes(x=poly_combined$logFC, y=-log10(poly_combined$P.Value), color = as.factor(to_color), label=poly_combined$external_gene_name)) + 
  facet_wrap(~ analysis) +
  geom_point(alpha=1, size = 2) +  scale_color_manual(labels = c('1' = paste0('FDR < 0.1'), 
                                                                 '2' = paste0('P-val < 0.005'), 
                                                                 '3' = paste0('n.s.')),
                                                      values=c('1' = "red", '2' = "orange", '3' = "grey60"))  + 
  guides(fill=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) + 
  ggtitle(paste0('Differential Gene Expression')) + xlab(expression('Effect Size (log'[2]*'(FC))')) + ylab(expression('\nSignificance (-log'[10]*'(P-Value))')) + labs(color = '') +
  geom_text_repel(data = top_poly_combined, aes(x=top_poly_combined$logFC, y=-log10(top_poly_combined$P.Value), label = top_poly_combined$external_gene_name), inherit.aes = FALSE,
                  size = 3, color="black", alpha=1, box.padding = unit(0.65, "lines")) +
  theme_bw() + theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
                     axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                     legend.text=element_text(size = 8),
                     legend.position = 'right',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     strip.text.x = element_text(size = 10),
                     legend.spacing.y = unit(1, "mm"),
                     legend.box.margin = margin(-10, -5, -10, -10))

# plot(poly_volcano_plot)

save(combined_volcano_plot, combined_GO_enrichment_plot, combined_EWCE_enrichment_plot, 
     poly_volcano_plot,
     file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')



### (5) Global logFC correlation of top poly1 vs poly2 DE genes ########################################################

# order both topTables by row name so all genes are in the same corresponding rows
poly1_ALL_ordered <- poly1_ALL[order(rownames(poly1_ALL)),]
poly2_ALL_ordered <- poly2_ALL[order(rownames(poly2_ALL)),]

# get the genes that are DE in either poly1 or poly2
poly1_ALL_DE <- poly1_ALL_ordered[which(poly1_ALL_ordered$P.Value < 0.005 | poly2_ALL_ordered$P.Value < 0.005),]
poly2_ALL_DE <- poly2_ALL_ordered[which(poly1_ALL_ordered$P.Value < 0.005 | poly2_ALL_ordered$P.Value < 0.005),]

# replace genes missing human homologs with their gene IDs
poly1_ALL_DE$external_gene_name <- ifelse(poly1_ALL_DE$external_gene_name == '', rownames(poly1_ALL_DE), poly1_ALL_DE$external_gene_name)
poly2_ALL_DE$external_gene_name <- ifelse(poly2_ALL_DE$external_gene_name == '', rownames(poly2_ALL_DE), poly2_ALL_DE$external_gene_name)

# create combined logFC data table and get colors for the data pronts
logFC <- poly1_ALL_DE[,c(8, 1)]
colnames(logFC) <- c('external_gene_name', 'poly1_logFC')
logFC['poly2_logFC'] <- poly2_ALL_DE$logFC
logFC['DGE'] <- ifelse(poly1_ALL_DE$P.Value < 0.005 & poly2_ALL_DE$P.Value < 0.005, 'Both', 
                       ifelse(poly1_ALL_DE$P.Value < 0.005, '1st Trim MIA', '2nd Trim MIA'))

# get the correlation between the logFC of DE genes
cor <- cor(logFC$poly1_logFC, logFC$poly2_logFC)

poly_logFC_plot <- ggplot(logFC, aes(x=logFC$poly1_logFC, y=logFC$poly2_logFC, color = logFC$DGE)) + 
  geom_point(size = 2) +
  xlim(-2, 2) + ylim(-2, 2) + geom_vline(xintercept=c(0), linetype="dashed") + geom_hline(yintercept=c(0), linetype="dashed") +
  geom_abline(intercept = 0, slope = -1, linetype = 'dashed') + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  annotate(geom = 'text', x = -1.25, y = 1.75, label = paste0("italic(R) ^ 2 == ", round(cor**2, 3)), parse = TRUE, size = 5) +
  ggtitle('Correlation of DE Genes') + 
  xlab(expression('1st Trim MIA (log'[2]*'(FC))')) + 
  ylab(expression('2nd Trim MIA (log'[2]*'(FC))')) +
  labs(color = 'P < 0.005 in:') +
  theme_bw() +
  geom_text_repel(
    data = subset(logFC, logFC$DGE == 'Both'),
    aes(x = subset(logFC, logFC$DGE == 'Both')$poly1_logFC, y = subset(logFC, logFC$DGE == 'Both')$poly2_logFC, label = subset(logFC, logFC$DGE == 'Both')$external_gene_name),
    size = 3,
    color = 'black',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = 'bottom',
        plot.margin = unit(c(5, 5, 5, 5), 'mm'),
        legend.margin=margin(c(0, 2, -5, 2), unit='mm'))
  
# plot(poly_logFC_plot)

save(combined_volcano_plot, combined_GO_enrichment_plot, combined_EWCE_enrichment_plot, 
     poly_volcano_plot, poly_logFC_plot,
     file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')



### (6) g:Profiler enrichment for top poly1 vs poly2 DE genes ##########################################################

library(gProfileR)

# order both topTables by row name so all genes are in the same corresponding rows
poly1_ALL_ordered <- poly1_ALL[order(rownames(poly1_ALL)),]
poly2_ALL_ordered <- poly2_ALL[order(rownames(poly2_ALL)),]

head(poly1_ALL_ordered)
head(poly2_ALL_ordered)

# get genes that are DE in poly1 but not changing at all in poly2
poly1_only <- poly1_ALL_ordered[which(poly1_ALL_ordered$P.Value < 0.05 & abs(poly2_ALL_ordered$logFC) < 0.2),]
poly1_only$signP <- -log10(poly1_only$P.Value) * sign(poly1_only$logFC)
poly1_only <- poly1_only[order(poly1_only$signP, decreasing = TRUE),]

poly1_only_UP <- rownames(poly1_only)[which(poly1_only$P.Value < 0.05 & poly1_only$logFC > 0)]
poly1_only_DOWN <- rev(rownames(poly1_only)[which(poly1_only$P.Value < 0.05 & poly1_only$logFC < 0)])


# get genes that are DE in poly2 but not changing at all in poly1
poly2_only <- poly2_ALL_ordered[which(poly2_ALL_ordered$P.Value < 0.05 & abs(poly1_ALL_ordered$logFC) < 0.2),]
poly2_only$signP <- -log10(poly2_only$P.Value) * sign(poly2_only$logFC)
poly2_only <- poly2_only[order(poly2_only$signP, decreasing = TRUE),]

poly2_only_UP <- rownames(poly2_only)[which(poly2_only$P.Value < 0.05 & poly2_only$logFC > 0)]
poly2_only_DOWN <- rev(rownames(poly2_only)[which(poly2_only$P.Value < 0.05 & poly2_only$logFC < 0)])


query_background <- rownames(mmul_CON_vs_MIA_ALL)

query_list <- list(poly1_only_UP, poly1_only_DOWN, poly2_only_UP, poly2_only_DOWN)
query_analysis_list <- list('poly1', 'poly1', 'poly2', 'poly2')
query_direction_list <- list('UP', 'DOWN', 'UP', 'DOWN')

GO_summary <- list()
top_GO_summary <- list()

mapply(function(query, query_analysis, query_direction) {
  
  # Perform g:Profiler enrichment based off of an ordered query for mmulatta, max set size = 1000
  go = gprofiler(query, organism="mmulatta", custom_bg = query_background,
                 correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
                 region_query = F, max_p_value = 1, min_set_size = 0, max_set_size=1000, numeric_ns = "",
                 include_graph = F,src_filter = c("GO", "KEGG"))
  GO_summary <<- rbind(GO_summary, data.frame('term' = go$term.name, 'p' = go$p.value, 'Analysis' = query_analysis, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  go = go[order(go$p.value)[1:min(4,nrow(go))],]
  top_GO_summary <<- rbind(top_GO_summary, data.frame('term' = go$term.name, 'p' = go$p.value, 'Analysis' = query_analysis, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  
  
}, query_list, query_analysis_list, query_direction_list)

GO_summary <- data.frame(GO_summary)
head(GO_summary)
hist(GO_summary$p)
GO_summary <- GO_summary[order(GO_summary$p, decreasing = FALSE),]
head(GO_summary)

top_GO_summary <- data.frame(top_GO_summary)
top_GO_summary

GO_summary <- top_GO_summary

library(ggplot2)
library(lemon)

GO_poly1_UP <- GO_summary[which(GO_summary$Direction == 'UP' & GO_summary$Analysis == 'poly1'),]
GO_poly1_UP$term <- c('nucleic acid phosphodiester\nbond hydrolysis',
                      'ncRNA processing',
                      'nucleosome',
                      'Elongator holoenzyme complex')

p1 <- ggplot(data = GO_poly1_UP, aes(x = GO_poly1_UP$term, y = -log10(GO_poly1_UP$p), fill = GO_poly1_UP$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_poly1_UP$term)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('p<0.05 in 1st Trim MIA & logFC<0.2 in 2nd Trim MIA') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_blank(),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10), 
        plot.margin = unit(c(2, 15, 2, 10),"mm"))

# plot(p1)


GO_poly1_DOWN <- GO_summary[which(GO_summary$Direction == 'DOWN' & GO_summary$Analysis == 'poly1'),]
GO_poly1_DOWN$term <- c('negative regulation of meiotic\nnuclear division',
                        'L-tyrosine aminotransferase activity',
                        'vitamin binding',
                        'carboxy-terminal domain protein\nkinase complex')

p2 <- ggplot(data = GO_poly1_DOWN, aes(x = GO_poly1_DOWN$term, y = -log10(GO_poly1_DOWN$p), fill = GO_poly1_DOWN$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_poly1_DOWN$term)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ylab(expression('\n-log'[10]*'(FDR)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10),
        plot.margin = unit(c(-2, 15, 2, 10),"mm"))

# plot(p2)


GO_poly2_UP <- GO_summary[which(GO_summary$Direction == 'UP' & GO_summary$Analysis == 'poly2'),]
GO_poly2_UP$term <- c('cellular response to cAMP',
                      'cell adhesion molecule binding',
                      'apical junction complex',
                      'heparin binding')

p3 <- ggplot(data = GO_poly2_UP, aes(x = GO_poly2_UP$term, y = -log10(GO_poly2_UP$p), fill = GO_poly2_UP$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_poly2_UP$term)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('p<0.05 in 2nd Trim MIA & logFC<0.2 in 1st Trim MIA') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_blank(),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10),
        plot.margin = unit(c(2, 15, 2, 10),"mm"))

# plot(p3)


GO_poly2_DOWN <- GO_summary[which(GO_summary$Direction == 'DOWN' & GO_summary$Analysis == 'poly2'),]
GO_poly2_DOWN$term <- c('RNA processing',
                        '90S preribosome',
                        'snoRNA binding',
                        'Y-form DNA binding')

p4 <- ggplot(data = GO_poly2_DOWN, aes(x = GO_poly2_DOWN$term, y = -log10(GO_poly2_DOWN$p), fill = GO_poly2_DOWN$Region)) +
  geom_bar(stat = 'identity', position = position_dodge(), fill = 'steelblue') +
  scale_x_discrete(limits = rev(GO_poly2_DOWN$term)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 6)) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ylab(expression('\n-log'[10]*'(FDR)')) +
  theme_bw() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text = element_blank(),
        legend.position = 'none',
        strip.text =  element_text(size = 10),
        plot.margin = unit(c(-2, 15, 2, 10),"mm"))

# plot(p4)


library(cowplot)

poly_GO_enrichment_plot <- plot_grid(p1, p2, p3, p4, ncol = 1, align = 'v')

# plot(poly_GO_enrichment_plot)

save(combined_volcano_plot, combined_GO_enrichment_plot, combined_EWCE_enrichment_plot, 
     poly_volcano_plot, poly_logFC_plot, poly_GO_enrichment_plot,
     file = './working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')




