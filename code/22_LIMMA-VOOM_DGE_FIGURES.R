#####################################################################################################################################
#
#   22_CONTE_NHP_MIA_DGE_FIGURES.R 
#
#   Generate figures for differential gene expression on 4 y/o non-human 
#   primate RNA-seq samples from DLPFC, ACC, HC, and V1 (Part 2)
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
library(edgeR)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

load('./working_data/filtered_unnormalized_datExpr_and_datMeta.RData')
load('./working_data/DGE_result_topTables.RData')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Volcano plot of regional differential gene expression
#   (2) Plot the total number of DE genes by region
#   (3) Cell-type enrichment from PsychEncode for global up and down genes with p < 0.005
#   (4) g:Profiler enrichment for top regional DE genes
#   (5) Heatmaps of cell-type specificity for top DE genes within each region
#   (6) Selected boxplots of top DE genes from each region
#   (7) Boxplots for upregulated myelin related genes in the hippocampus
#   (8) logFC correlation for all genes across 4 regions
#   (9) logFC correlation for genes with p < 0.005 across 4 regions
#
########################################################################################################################

### (1) Volcano plot of regional differential gene expression ##########################################################

mmul_CON_vs_MIA_DLPFC <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
mmul_CON_vs_MIA_DLPFC['Region'] <- 1

mmul_CON_vs_MIA_ACC <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
mmul_CON_vs_MIA_ACC['Region'] <- 2

mmul_CON_vs_MIA_HC <- DGE_result_topTables$mmul_CON_vs_MIA_HC
mmul_CON_vs_MIA_HC['Region'] <- 3

mmul_CON_vs_MIA_V1 <- DGE_result_topTables$mmul_CON_vs_MIA_V1
mmul_CON_vs_MIA_V1['Region'] <- 4

mmul_CON_vs_MIA_Combined <- rbind(mmul_CON_vs_MIA_DLPFC,
                                  mmul_CON_vs_MIA_ACC,
                                  mmul_CON_vs_MIA_HC,
                                  mmul_CON_vs_MIA_V1)


mmul_CON_vs_MIA_Combined['all_FDR'] <- fdrtool::fdrtool(mmul_CON_vs_MIA_Combined$t, plot = TRUE)$qval

top_mmul_CON_vs_MIA_Combined <- rbind(mmul_CON_vs_MIA_DLPFC[1:10,],
                                      mmul_CON_vs_MIA_ACC[1:10,],
                                      mmul_CON_vs_MIA_HC[1:10,],
                                      mmul_CON_vs_MIA_V1[1:10,])

# use and ifelse statement to assign the correct color to genes in the volcano plot based on their P-Values
to_color <- c(ifelse(mmul_CON_vs_MIA_Combined$all_FDR < 0.1, 1, ifelse(mmul_CON_vs_MIA_Combined$P.Value < 0.005 & mmul_CON_vs_MIA_Combined$all_FDR >= 0.1 ,2, 3)))
to_color_table <- table(to_color)

region.labels <- c('DLPFC', 'ACC', 'HC', 'V1')
names(region.labels) <- c(1, 2, 3, 4)

regional_volcano_plot <- ggplot(as.data.frame(mmul_CON_vs_MIA_Combined), aes(x=mmul_CON_vs_MIA_Combined$logFC, y=-log10(mmul_CON_vs_MIA_Combined$P.Value), color = as.factor(to_color), label=mmul_CON_vs_MIA_Combined$external_gene_name)) + 
  facet_wrap(~ Region, nrow = 2, labeller = labeller(Region = region.labels)) +
  geom_point(alpha=1, size = 2) +  scale_color_manual(labels = c('1' = paste0('FDR < 0.1'), 
                                                                 '2' = paste0('P-val < 0.005'), 
                                                                 '3' = paste0('n.s.')),
                                                      values=c('1' = "red", '2' = "orange", '3' = "grey60"))  + 
  guides(fill=guide_legend(keywidth=0.1,keyheight=0.5,default.unit="inch")) + 
  ggtitle(paste0('Regional Differential Gene Expression')) + xlab(expression('Effect Size (log'[2]*'(FC))')) + ylab(expression('\nSignificance (-log'[10]*'(P-Value))')) + labs(color = '') +
  geom_text_repel(data = top_mmul_CON_vs_MIA_Combined, aes(x=top_mmul_CON_vs_MIA_Combined$logFC, y=-log10(top_mmul_CON_vs_MIA_Combined$P.Value), label = top_mmul_CON_vs_MIA_Combined$external_gene_name), inherit.aes = FALSE,
                  size = 3, color="black", alpha=1, box.padding = unit(0.65, "lines")) +
  theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                     axis.title.x = element_text(size = 10, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                     legend.text = element_text(size = 8),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 10),
                     strip.text.x = element_text(size = 10),
                     legend.margin = margin(0, 0, 0, 0),
                     legend.box.margin = margin(-10, -10, -10, -10),
                     legend.position = 'bottom',
                     legend.direction = 'horizontal',
                     plot.margin = unit(c(2, 2, 2, 2), 'mm'))

# plot(regional_volcano_plot)

save(regional_volcano_plot, file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (2) Plot the total number of DE genes by region ####################################################################

library(plyr)
mmul_CON_vs_MIA_Combined$Region <- mapvalues(mmul_CON_vs_MIA_Combined$Region, from = c(1, 2, 3, 4), to = c('DLPFC', 'ACC', 'HC', 'V1'))

DLPFC_005 <- length(which(mmul_CON_vs_MIA_Combined$Region == 'DLPFC' & mmul_CON_vs_MIA_Combined$P.Value < 0.005))
ACC_005 <- length(which(mmul_CON_vs_MIA_Combined$Region == 'ACC' & mmul_CON_vs_MIA_Combined$P.Value < 0.005))
HC_005 <- length(which(mmul_CON_vs_MIA_Combined$Region == 'HC' & mmul_CON_vs_MIA_Combined$P.Value < 0.005))
V1_005 <- length(which(mmul_CON_vs_MIA_Combined$Region == 'V1' & mmul_CON_vs_MIA_Combined$P.Value < 0.005))

DLPFC_FDR <- length(which(mmul_CON_vs_MIA_Combined$Region == 'DLPFC' & mmul_CON_vs_MIA_Combined$all_FDR < 0.1))
ACC_FDR <- length(which(mmul_CON_vs_MIA_Combined$Region == 'ACC' & mmul_CON_vs_MIA_Combined$all_FDR < 0.1))
HC_FDR <- length(which(mmul_CON_vs_MIA_Combined$Region == 'HC' & mmul_CON_vs_MIA_Combined$all_FDR < 0.1))
V1_FDR <- length(which(mmul_CON_vs_MIA_Combined$Region == 'V1' & mmul_CON_vs_MIA_Combined$all_FDR < 0.1))

num_DE_genes_005 <- c(DLPFC_005, ACC_005, HC_005, V1_005)
names(num_DE_genes_005) <- c('DLPFC', 'ACC', 'HC', 'V1')

num_DE_genes_FDR <- c(DLPFC_FDR, ACC_FDR, HC_FDR, V1_FDR)
names(num_DE_genes_FDR) <- c('DLPFC', 'ACC', 'HC', 'V1')

num_DE_genes <- data.frame('num' = c(num_DE_genes_005, num_DE_genes_FDR),
                           'sig' = c(rep('p<0.005', 4), rep('FDR<0.1', 4)),
                           'region' = rep(c('DLPFC', 'ACC', 'HC', 'V1'), 2))

num_DE_genes$sig <- factor(num_DE_genes$sig)
num_DE_genes$region <- factor(num_DE_genes$region)

num_DE_genes$sig <- factor(num_DE_genes$sig, levels = levels(num_DE_genes$sig)[c(2,1)])
num_DE_genes$region <- factor(num_DE_genes$region, levels = levels(num_DE_genes$region)[c(2,1,3,4)])

num_DE_genes_plot <- ggplot(data = num_DE_genes, aes(x = num_DE_genes$region, y = num_DE_genes$num, fill = num_DE_genes$sig)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label=num_DE_genes$num), vjust = -0.2, color="black",
            position = position_dodge(0.9), size = 3) +
  scale_fill_manual(values = c('orange','red')) +
  ggtitle('# DE Genes by Region') +
  ylab('# DE Genes\n') +
  ylim(c(0, 325)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(5, "mm"),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(2, 2, -5, 2), 'mm'))

# plot(num_DE_genes_plot)

save(regional_volcano_plot, num_DE_genes_plot, file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (3) Cell-type enrichment from PsychEncode for global up and down genes with p < 0.005 ##############################

library(EWCE)

load('./working_data/CellTypeData_PEC_Lake_Adult_Human_nucSeq.rda')

dim(ctd[[2]]$specificity)
head(ctd[[2]]$specificity)

load('./working_data/DGE_result_topTables.RData')

mmul_CON_vs_MIA_ALL <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
mmul_CON_vs_MIA_DLPFC <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
mmul_CON_vs_MIA_ACC <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
mmul_CON_vs_MIA_HC <- DGE_result_topTables$mmul_CON_vs_MIA_HC
mmul_CON_vs_MIA_V1 <- DGE_result_topTables$mmul_CON_vs_MIA_V1

DLPFC_UP = mmul_CON_vs_MIA_DLPFC$external_gene_name[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.05 & mmul_CON_vs_MIA_DLPFC$logFC > 0)]
DLPFC_DOWN = mmul_CON_vs_MIA_DLPFC$external_gene_name[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.05 & mmul_CON_vs_MIA_DLPFC$logFC < 0)]
ACC_UP = mmul_CON_vs_MIA_ACC$external_gene_name[which(mmul_CON_vs_MIA_ACC$P.Value < 0.05 & mmul_CON_vs_MIA_ACC$logFC > 0)]
ACC_DOWN = mmul_CON_vs_MIA_ACC$external_gene_name[which(mmul_CON_vs_MIA_ACC$P.Value < 0.05 & mmul_CON_vs_MIA_ACC$logFC < 0)]
HC_UP = mmul_CON_vs_MIA_HC$external_gene_name[which(mmul_CON_vs_MIA_HC$P.Value < 0.05 & mmul_CON_vs_MIA_HC$logFC > 0)]
HC_DOWN = mmul_CON_vs_MIA_HC$external_gene_name[which(mmul_CON_vs_MIA_HC$P.Value < 0.05 & mmul_CON_vs_MIA_HC$logFC < 0)]
V1_UP = mmul_CON_vs_MIA_V1$external_gene_name[which(mmul_CON_vs_MIA_V1$P.Value < 0.05 & mmul_CON_vs_MIA_V1$logFC > 0)]
V1_DOWN = mmul_CON_vs_MIA_V1$external_gene_name[which(mmul_CON_vs_MIA_V1$P.Value < 0.05 & mmul_CON_vs_MIA_V1$logFC < 0)]


set.seed(12211991)
reps = 10000

background_adult = rownames(ctd[[2]]$specificity)
length(background_adult)
background_adult = background_adult[which(background_adult %in% mmul_CON_vs_MIA_ALL$external_gene_name)]
length(background_adult)

DLPFC_UP = DLPFC_UP[which(DLPFC_UP %in% background_adult)]
DLPFC_DOWN = DLPFC_DOWN[which(DLPFC_DOWN %in% background_adult)]
ACC_UP = ACC_UP[which(ACC_UP %in% background_adult)]
ACC_DOWN = ACC_DOWN[which(ACC_DOWN %in% background_adult)]
HC_UP = HC_UP[which(HC_UP %in% background_adult)]
HC_DOWN = HC_DOWN[which(HC_DOWN %in% background_adult)]
V1_UP = V1_UP[which(V1_UP %in% background_adult)]
V1_DOWN = V1_DOWN[which(V1_DOWN %in% background_adult)]


# Bootstrap significance testing, without controlling for transcript length and GC content
DLPFC_UP_results = bootstrap.enrichment.test(sct_data=ctd, hits=DLPFC_UP, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
DLPFC_DOWN_results = bootstrap.enrichment.test(sct_data=ctd, hits=DLPFC_DOWN, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
ACC_UP_results = bootstrap.enrichment.test(sct_data=ctd, hits=ACC_UP, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
ACC_DOWN_results = bootstrap.enrichment.test(sct_data=ctd, hits=ACC_DOWN, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
HC_UP_results = bootstrap.enrichment.test(sct_data=ctd, hits=HC_UP, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
HC_DOWN_results = bootstrap.enrichment.test(sct_data=ctd, hits=HC_DOWN, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
V1_UP_results = bootstrap.enrichment.test(sct_data=ctd, hits=V1_UP, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")
V1_DOWN_results = bootstrap.enrichment.test(sct_data=ctd, hits=V1_DOWN, bg=background_adult,reps=reps,annotLevel=2,genelistSpecies="human",sctSpecies="human")

ewce_results <- rbind(data.frame(DLPFC_UP_results$results[,c(1,3,5)], 'Region' = '1', 'Direction' = '1'),
                      data.frame(DLPFC_DOWN_results$results[,c(1,3,5)], 'Region' = '1', 'Direction' = '2'),
                      data.frame(ACC_UP_results$results[,c(1,3,5)], 'Region' = '2', 'Direction' = '1'),
                      data.frame(ACC_DOWN_results$results[,c(1,3,5)], 'Region' = '2', 'Direction' = '2'),
                      data.frame(HC_UP_results$results[,c(1,3,5)], 'Region' = '3', 'Direction' = '1'),
                      data.frame(HC_DOWN_results$results[,c(1,3,5)], 'Region' = '3', 'Direction' = '2'),
                      data.frame(V1_UP_results$results[,c(1,3,5)], 'Region' = '4', 'Direction' = '1'),
                      data.frame(V1_DOWN_results$results[,c(1,3,5)], 'Region' = '4', 'Direction' = '2'))

ewce_results$FDR <- p.adjust(ewce_results$p, method = 'fdr')
ewce_results$sd_from_mean <- ifelse(ewce_results$sd_from_mean > 0, ewce_results$sd_from_mean, 0)
print(ewce_results)

ewce_results$Direction <- factor(ewce_results$Direction, labels = c('UP', 'DOWN'))
ewce_results$Region <- factor(ewce_results$Region, labels = c('DLPFC', 'ACC', 'HC', 'V1'))
ewce_results$sig <- ifelse(ewce_results$FDR < 0.05, '*', '')

print(ewce_results)

regional_EWCE_enrichment_plot <- ggplot(data = ewce_results, aes(x = ewce_results$CellType, y = ewce_results$sd_from_mean, fill = ewce_results$Region)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_grid(rows = vars(Direction), cols = vars(Region)) +
  geom_text(aes(label=ewce_results$sig), vjust=0.5, color="red",
            position = position_dodge(0), size = 6) +
  ggtitle('Cell-type Enrichment from PsychEncode') +
  ylab('Standard Deviation from Mean') +
  ylim(c(0, 25)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.key.size = unit(5, "mm"),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 8),
        strip.text =  element_text(size = 10),
        plot.margin=unit(c(10,2,5,2),"mm"))


plot(regional_EWCE_enrichment_plot)

save(regional_volcano_plot, num_DE_genes_plot, file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (4) g:Profiler enrichment for top regional DE genes ################################################################

query_background <- rownames(mmul_CON_vs_MIA_ALL)

# get top up and down gene lists for each region
DLPFC_UP = rownames(mmul_CON_vs_MIA_DLPFC)[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.05 & mmul_CON_vs_MIA_DLPFC$logFC > 0)]
DLPFC_DOWN = rev(rownames(mmul_CON_vs_MIA_DLPFC)[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.05 & mmul_CON_vs_MIA_DLPFC$logFC < 0)])
ACC_UP = rownames(mmul_CON_vs_MIA_ACC)[which(mmul_CON_vs_MIA_ACC$P.Value < 0.05 & mmul_CON_vs_MIA_ACC$logFC > 0)]
ACC_DOWN = rev(rownames(mmul_CON_vs_MIA_ACC)[which(mmul_CON_vs_MIA_ACC$P.Value < 0.05 & mmul_CON_vs_MIA_ACC$logFC < 0)])
HC_UP = rownames(mmul_CON_vs_MIA_HC)[which(mmul_CON_vs_MIA_HC$P.Value < 0.05 & mmul_CON_vs_MIA_HC$logFC > 0)]
HC_DOWN = rev(rownames(mmul_CON_vs_MIA_HC)[which(mmul_CON_vs_MIA_HC$P.Value < 0.05 & mmul_CON_vs_MIA_HC$logFC < 0)])
V1_UP = rownames(mmul_CON_vs_MIA_V1)[which(mmul_CON_vs_MIA_V1$P.Value < 0.05 & mmul_CON_vs_MIA_V1$logFC > 0)]
V1_DOWN = rev(rownames(mmul_CON_vs_MIA_V1)[which(mmul_CON_vs_MIA_V1$P.Value < 0.05 & mmul_CON_vs_MIA_V1$logFC < 0)])

# remove genes lacking a human homolog
DLPFC_UP <- DLPFC_UP[which(DLPFC_UP != '')]
DLPFC_DOWN <- DLPFC_DOWN[which(DLPFC_DOWN != '')]
ACC_UP <- ACC_UP[which(ACC_UP != '')]
ACC_DOWN <- ACC_DOWN[which(ACC_DOWN != '')]
HC_UP <- HC_UP[which(HC_UP != '')]
HC_DOWN <- HC_DOWN[which(HC_DOWN != '')]
V1_UP <- V1_UP[which(V1_UP != '')]
V1_DOWN <- V1_DOWN[which(V1_DOWN != '')]

query_list <- list(DLPFC_UP, DLPFC_DOWN, ACC_UP, ACC_DOWN, HC_UP, HC_DOWN, V1_UP, V1_DOWN)
query_region_list <- list('DLPFC', 'DLPFC', 'ACC', 'ACC', 'HC', 'HC', 'V1', 'V1')
query_direction_list <- list('UP', 'DOWN', 'UP', 'DOWN', 'UP', 'DOWN', 'UP', 'DOWN')

GO_summary <- list()
top_GO_summary <- list()

mapply(function(query, query_region, query_direction) {
  
  go = gprofiler(query, organism="mmulatta", custom_bg = query_background,
                 correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
                 region_query = F, max_p_value = 1, min_set_size = 0, max_set_size = 1000, numeric_ns = "",
                 include_graph = F, src_filter = c("GO", "KEGG"))
  GO_summary <<- rbind(GO_summary, data.frame('term' = go$term.name, 'p' = go$p.value, 'Region' = query_region, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  go = go[order(go$p.value)[1:min(3,nrow(go))],]
  top_GO_summary <<- rbind(top_GO_summary, data.frame('term' = go$term.name, 'p' = go$p.value, 'Region' = query_region, 'Direction' = query_direction, 'core_enrichment' = go$intersection))
  
  
}, query_list, query_region_list, query_direction_list)

GO_summary <- data.frame(GO_summary)
head(GO_summary)
hist(GO_summary$p)


top_GO_summary <- data.frame(top_GO_summary)
top_GO_summary

GO_summary <- top_GO_summary


library(ggplot2)
library(lemon)

GO_summary_UP <- GO_summary[which(GO_summary$Direction == 'UP'),]
GO_summary_UP$term <- c('protein kinase B signaling',
                        'protein tyrosine kinase activity',
                        'protein membrane anchor',
                        'cellular response to UV-C',
                        'urokinase plasminogen activator\nsignaling pathway',
                        'maintenance of epithelial cell\napical/basal polarity',
                        'axon ensheathment',
                        'extracellular space',
                        'transporter activity',
                        'mitochondrial part',
                        'protein repair',
                        'peptidyl-lysine deglutarylation')

GO_summary_UP$term <- factor(GO_summary_UP$term, levels = rev(GO_summary_UP$term[c(1:12)]))
GO_summary_UP$Region <- factor(GO_summary_UP$Region, levels = c('DLPFC', 'ACC', 'HC', 'V1'))

p1 <- ggplot(data = GO_summary_UP, aes(x = GO_summary_UP$term, y = -log10(GO_summary_UP$p), fill = GO_summary_UP$Region)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  ylim(c(0, 16)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_blank(),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 8),
        strip.text =  element_text(size = 10),
        plot.margin = unit(c(10, 2, 2, 2),"mm"))

# plot(p1)


GO_summary_DOWN <- GO_summary[which(GO_summary$Direction == 'DOWN'),]
GO_summary_DOWN$term <- c('structural constituent of\nribosome',
                          'ribosome',
                          'translation',
                          'positive regulation of\nacrosome reaction',
                          'Hsp90 protein binding',
                          'IRE1-RACK1-PP2A complex',
                          'anterograde trans-synaptic\nsignaling',
                          'synapse',
                          'ion gated channel activity',
                          'chromosome organization',
                          'nucleoplasm part',
                          'chromatin binding')

# GO_summary_DOWN$term <- factor(GO_summary_DOWN$term, levels = rev(GO_summary_DOWN$term[c(3,2,1,6,5,4,9,8,7,12,11,10)]))
GO_summary_DOWN$term <- factor(GO_summary_DOWN$term, levels = rev(GO_summary_DOWN$term))
GO_summary_DOWN$Region <- factor(GO_summary_DOWN$Region, levels = c('DLPFC', 'ACC', 'HC', 'V1'))

p2 <- ggplot(data = GO_summary_DOWN, aes(x = GO_summary_DOWN$term, y = -log10(GO_summary_DOWN$p), fill = GO_summary_DOWN$Region)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_rep_wrap(~ Direction, nrow = 1, strip.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  ylim(c(0, 16)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_text(size = 10, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 8),
        strip.text =  element_text(size = 10),
        legend.margin = margin(0, 0, 0, -20),
        legend.box.margin = margin(-2, 0, 0, -20),
        plot.margin=unit(c(-2,5,2,5),"mm"))

# plot(p2)


library(cowplot)

#pdf('./figures/Regional GO Enrichment.pdf', 8, 6)
regional_GO_enrichment_plot <- plot_grid(p1, p2, ncol = 1, align = 'v')
#dev.off()

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (5) Heatmaps of cell-type specificity for top DE genes within each region ##########################################

library(reshape2)

head(ctd[[2]]$specificity)

DLPFC_UP = mmul_CON_vs_MIA_DLPFC$external_gene_name[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.005 & mmul_CON_vs_MIA_DLPFC$logFC > 0)]
DLPFC_DOWN = mmul_CON_vs_MIA_DLPFC$external_gene_name[which(mmul_CON_vs_MIA_DLPFC$P.Value < 0.005 & mmul_CON_vs_MIA_DLPFC$logFC < 0)]
ACC_UP = mmul_CON_vs_MIA_ACC$external_gene_name[which(mmul_CON_vs_MIA_ACC$P.Value < 0.005 & mmul_CON_vs_MIA_ACC$logFC > 0)]
ACC_DOWN = mmul_CON_vs_MIA_ACC$external_gene_name[which(mmul_CON_vs_MIA_ACC$P.Value < 0.005 & mmul_CON_vs_MIA_ACC$logFC < 0)]
HC_UP = mmul_CON_vs_MIA_HC$external_gene_name[which(mmul_CON_vs_MIA_HC$P.Value < 0.005 & mmul_CON_vs_MIA_HC$logFC > 0)]
HC_DOWN = mmul_CON_vs_MIA_HC$external_gene_name[which(mmul_CON_vs_MIA_HC$P.Value < 0.005 & mmul_CON_vs_MIA_HC$logFC < 0)]
V1_UP = mmul_CON_vs_MIA_V1$external_gene_name[which(mmul_CON_vs_MIA_V1$P.Value < 0.005 & mmul_CON_vs_MIA_V1$logFC > 0)]
V1_DOWN = mmul_CON_vs_MIA_V1$external_gene_name[which(mmul_CON_vs_MIA_V1$P.Value < 0.005 & mmul_CON_vs_MIA_V1$logFC < 0)]


DLPFC_UP_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% DLPFC_UP),]
DLPFC_UP_specificity <- DLPFC_UP_specificity[order(apply(DLPFC_UP_specificity, 1, max), decreasing = TRUE),]
DLPFC_UP_specificity <- DLPFC_UP_specificity[1:10,]

DLPFC_UP_specificity <- melt(DLPFC_UP_specificity)
colnames(DLPFC_UP_specificity) <- c('gene', 'cell', 'specificity')
DLPFC_UP_specificity$specificity <- DLPFC_UP_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p1 <- ggplot(DLPFC_UP_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(DLPFC_UP_specificity$specificity > 50, round(DLPFC_UP_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  ggtitle('DLPFC') + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin=unit(c(2,0,-10,0),"mm"))

plot(p1)


###


DLPFC_DOWN_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% DLPFC_DOWN),]
DLPFC_DOWN_specificity <- DLPFC_DOWN_specificity[order(apply(DLPFC_DOWN_specificity, 1, max), decreasing = TRUE),]
DLPFC_DOWN_specificity <- DLPFC_DOWN_specificity[1:10,]

DLPFC_DOWN_specificity <- melt(DLPFC_DOWN_specificity)
colnames(DLPFC_DOWN_specificity) <- c('gene', 'cell', 'specificity')
DLPFC_DOWN_specificity$specificity <- DLPFC_DOWN_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p2 <- ggplot(DLPFC_DOWN_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(DLPFC_DOWN_specificity$specificity > 50, round(DLPFC_DOWN_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(1, "cm"),
        legend.position = 'none',
        legend.title = element_blank(),
        plot.margin=unit(c(0,0,5,0),"mm"))

plot(p2)

###

ACC_UP_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% ACC_UP),]
ACC_UP_specificity <- ACC_UP_specificity[order(apply(ACC_UP_specificity, 1, max), decreasing = TRUE),]
ACC_UP_specificity <- ACC_UP_specificity[1:4,]

ACC_UP_specificity <- melt(ACC_UP_specificity)
colnames(ACC_UP_specificity) <- c('gene', 'cell', 'specificity')
ACC_UP_specificity$specificity <- ACC_UP_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p3 <- ggplot(ACC_UP_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(ACC_UP_specificity$specificity > 50, round(ACC_UP_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  ggtitle('ACC') + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin=unit(c(2,0,-10,0),"mm"))

plot(p3)


###


ACC_DOWN_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% ACC_DOWN),]
ACC_DOWN_specificity <- ACC_DOWN_specificity[order(apply(ACC_DOWN_specificity, 1, max), decreasing = TRUE),]
ACC_DOWN_specificity <- ACC_DOWN_specificity[1:6,]

ACC_DOWN_specificity <- melt(ACC_DOWN_specificity)
colnames(ACC_DOWN_specificity) <- c('gene', 'cell', 'specificity')
ACC_DOWN_specificity$specificity <- ACC_DOWN_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p4 <- ggplot(ACC_DOWN_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(ACC_DOWN_specificity$specificity > 50, round(ACC_DOWN_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(1, "cm"),
        legend.position = 'none',
        legend.title = element_blank(),
        plot.margin=unit(c(0,0,5,0),"mm"))

plot(p4)


###


HC_UP_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% HC_UP),]
HC_UP_specificity <- HC_UP_specificity[order(apply(HC_UP_specificity, 1, max), decreasing = TRUE),]
HC_UP_specificity <- HC_UP_specificity[1:10,]

HC_UP_specificity <- melt(HC_UP_specificity)
colnames(HC_UP_specificity) <- c('gene', 'cell', 'specificity')
HC_UP_specificity$specificity <- HC_UP_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p5 <- ggplot(HC_UP_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(HC_UP_specificity$specificity > 50, round(HC_UP_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  ggtitle('HC') + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin=unit(c(2,0,-10,0),"mm"))

plot(p5)


###


HC_DOWN_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% HC_DOWN),]
HC_DOWN_specificity <- HC_DOWN_specificity[order(apply(HC_DOWN_specificity, 1, max), decreasing = TRUE),]
HC_DOWN_specificity <- HC_DOWN_specificity[1:10,]

HC_DOWN_specificity <- melt(HC_DOWN_specificity)
colnames(HC_DOWN_specificity) <- c('gene', 'cell', 'specificity')
HC_DOWN_specificity$specificity <- HC_DOWN_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p6 <- ggplot(HC_DOWN_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(HC_DOWN_specificity$specificity > 50, round(HC_DOWN_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(1, "cm"),
        legend.position = 'none',
        legend.title = element_blank(),
        plot.margin=unit(c(0,0,5,0),"mm"))

plot(p6)


###


V1_UP_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% V1_UP),]
V1_UP_specificity <- V1_UP_specificity[order(apply(V1_UP_specificity, 1, max), decreasing = TRUE),]
V1_UP_specificity <- V1_UP_specificity[1:10,]

V1_UP_specificity <- melt(V1_UP_specificity)
colnames(V1_UP_specificity) <- c('gene', 'cell', 'specificity')
V1_UP_specificity$specificity <- V1_UP_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p7 <- ggplot(V1_UP_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(V1_UP_specificity$specificity > 50, round(V1_UP_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  ggtitle('V1') + 
  theme_bw() + 
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 0, b = 2, l = 0)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin=unit(c(2,0,-10,0),"mm"))

plot(p7)


###

V1_DOWN_specificity <- ctd[[2]]$specificity[which(rownames(ctd[[2]]$specificity) %in% V1_DOWN),]
V1_DOWN_specificity <- V1_DOWN_specificity[order(apply(V1_DOWN_specificity, 1, max), decreasing = TRUE),]
V1_DOWN_specificity <- V1_DOWN_specificity[1:10,]

V1_DOWN_specificity <- melt(V1_DOWN_specificity)
colnames(V1_DOWN_specificity) <- c('gene', 'cell', 'specificity')
V1_DOWN_specificity$specificity <- V1_DOWN_specificity$specificity * 100


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

p8 <- ggplot(V1_DOWN_specificity, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(V1_DOWN_specificity$specificity > 50, round(V1_DOWN_specificity$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5, 0, 0, 0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(1, "cm"),
        legend.position = 'none',
        legend.direction = 'horizontal',
        legend.title = element_text(size = 8),
        plot.margin=unit(c(-10,0,5,0),"mm"))

plot(p8)


###

plot_list <- list(p1, p3, p5, p7, p2, p4, p6, p8)
cell_type_heatmaps <- plot_grid(plotlist = plot_list, align = 'hv', nrow = 2, ncol = 4)

###

library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

test_heatmap <- data.frame('cell' = c(1, 2), 'gene' = c(1, 2), 'specificity' = c(0, 100))

specificity_scale_bar <- ggplot(test_heatmap, aes(x = cell, y = gene, fill = specificity)) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(test_heatmap$specificity > 50, round(test_heatmap$specificity, 1), '')), 
            size = 1.6, color = 'white', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 100)) +
  theme_bw() + 
  labs(fill = '% Cell-Type Specificity') +
  #guides(fill = guide_legend(title.position = "top")) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5, 0, 0, 0)),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.margin=unit(c(-10,0,5,0),"mm"))

# plot(specificity_scale_bar)

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     cell_type_heatmaps, specificity_scale_bar,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (6) Selected boxplots of top DE genes from each region #############################################################

library(scales)

tt.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
DLPFC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_V1

datExpr <- mmul_txi.gene
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr)), method = 'TMM', log = TRUE)

genes <- c('RFC2', 'DGCR2', 'MLC1', 'KDSR', 'PACS1', 'PGBD3', 'GRM3', 'GRM8', 'TRIO', 'GRM1', 'GRIK2', 'GRID2')
gene_ids <- mmul_annotEns87$ensembl_gene_id[match(genes, mmul_annotEns87$hsapiens_homolog_associated_gene_name)]


RFC2 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[1]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[1]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[1], subtitle = paste0('DLPFC: p=', round(DLPFC.fullModel$P.Value[which(DLPFC.fullModel$external_gene_name == genes[1])], 4))) +
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

DGCR2 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[2]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[2]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[2], subtitle = paste0('DLPFC: p=', round(DLPFC.fullModel$P.Value[which(DLPFC.fullModel$external_gene_name == genes[2])], 4))) +
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

MLC1 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[3]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[3]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[3], subtitle = paste0('DLPFC: p=', round(DLPFC.fullModel$P.Value[which(DLPFC.fullModel$external_gene_name == genes[3])], 4))) +
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

KDSR <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[4]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[4]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[4], subtitle = paste0('ACC: p=', round(ACC.fullModel$P.Value[which(ACC.fullModel$external_gene_name == genes[4])], 4))) +
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

PACS1 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[5]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[5]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[5], subtitle = paste0('ACC: p=', round(ACC.fullModel$P.Value[which(ACC.fullModel$external_gene_name == genes[5])], 4))) +
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

PGBD3 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[6]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[6]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[6], subtitle = paste0('ACC: p=', round(ACC.fullModel$P.Value[which(ACC.fullModel$external_gene_name == genes[6])], 4))) +
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

GRM3 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[7]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[7]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[7], subtitle = paste0('HC: p=', round(HC.fullModel$P.Value[which(HC.fullModel$external_gene_name == genes[7])], 4))) +
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

GRM8 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[8]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[8]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[8], subtitle = paste0('HC: p=', round(HC.fullModel$P.Value[which(HC.fullModel$external_gene_name == genes[8])], 4))) +
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

TRIO <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[9]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[9]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[9], subtitle = paste0('HC: p=', round(HC.fullModel$P.Value[which(HC.fullModel$external_gene_name == genes[9])], 4))) +
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

GRM1 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[10]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[10]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[10], subtitle = paste0('V1: p=', round(V1.fullModel$P.Value[which(V1.fullModel$external_gene_name == genes[10])], 4))) +
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

GRIK2 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[11]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[11]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[11], subtitle = paste0('V1: p=', round(V1.fullModel$P.Value[which(V1.fullModel$external_gene_name == genes[11])], 4))) +
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

GRID2 <- ggplot(data.frame(datExpr.norm[which(rownames(datExpr.norm) == gene_ids[12]),]), aes(x=datMeta$MIA, y=datExpr.norm[which(rownames(datExpr.norm) == gene_ids[12]),], color = datMeta$MIA)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.6) +
  #geom_jitter(shape=16, position=position_jitter(0.1), size = 4) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = 0.75)) +
  facet_grid(. ~ datMeta$Region) + 
  ggtitle(genes[12], subtitle = paste0('V1: p=', round(V1.fullModel$P.Value[which(V1.fullModel$external_gene_name == genes[12])], 4))) +
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

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     cell_type_heatmaps, specificity_scale_bar, RFC2, DGCR2, MLC1, KDSR,
     PACS1, PGBD3, GRM3, GRM8, TRIO, GRM1, GRIK2, GRID2,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (7) Boxplots for upregulated mylein related genes in the hippocampus ###############################################

library(scales)

to_keep <- which(datMeta$Region == 'HC')
HC_datMeta <- datMeta[to_keep,]

mmul_HC.fullModel <- DGE_result_topTables$mmul_CON_vs_MIA_HC
mmul_datExpr <- mmul_txi.gene
mmul_datExpr.norm <- cpm(calcNormFactors(DGEList(mmul_datExpr)), method = 'TMM', log = TRUE)
mmul_datExpr.norm <- mmul_datExpr.norm[,to_keep] # has 14975 genes total

hg_HC.fullModel <- DGE_result_topTables$hg_CON_vs_MIA_HC
hg_datExpr <- hg_txi.gene
hg_datExpr.norm <- cpm(calcNormFactors(DGEList(hg_datExpr)), method = 'TMM', log = TRUE)
hg_datExpr.norm <- hg_datExpr.norm[,to_keep] # has 17640 genes total

myelin_genes <- c('OLIG2', 'SOX10', 'MYRF', 'PLLP', 'MBP', 'CNP', 'MOG', 'MAG')

OLIG2 <- data.frame(rep(NA, 13))
OLIG2$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'OLIG2')]),]
OLIG2[,1] <- OLIG2$expr/mean(OLIG2$expr[which(HC_datMeta$MIA == 'CON')])
colnames(OLIG2) <- c('norm', 'expr')

SOX10 <- data.frame(rep(NA, 13))
SOX10$expr <- hg_datExpr.norm[which(rownames(hg_datExpr.norm) == rownames(hg_HC.fullModel)[which(hg_HC.fullModel$external_gene_name == 'SOX10')]),]
SOX10[,1] <- SOX10$expr/mean(SOX10$expr[which(HC_datMeta$MIA == 'CON')])
colnames(SOX10) <- c('norm', 'expr')

MYRF <- data.frame(rep(NA, 13))
MYRF$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'MYRF')]),]
MYRF[,1] <- MYRF$expr/mean(MYRF$expr[which(HC_datMeta$MIA == 'CON')])
colnames(MYRF) <- c('norm', 'expr')

PLLP <- data.frame(rep(NA, 13))
PLLP$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'PLLP')]),]
PLLP[,1] <- PLLP$expr/mean(PLLP$expr[which(HC_datMeta$MIA == 'CON')])
colnames(PLLP) <- c('norm', 'expr')

MBP <- data.frame(rep(NA, 13))
MBP$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'MBP')]),]
MBP[,1] <- MBP$expr/mean(MBP$expr[which(HC_datMeta$MIA == 'CON')])
colnames(MBP) <- c('norm', 'expr')

CNP <- data.frame(rep(NA, 13))
CNP$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'CNP')]),]
CNP[,1] <- CNP$expr/mean(CNP$expr[which(HC_datMeta$MIA == 'CON')])
colnames(CNP) <- c('norm', 'expr')

MOG <- data.frame(rep(NA, 13))
MOG$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'MOG')]),]
MOG[,1] <- MOG$expr/mean(MOG$expr[which(HC_datMeta$MIA == 'CON')])
colnames(MOG) <- c('norm', 'expr')

MAG <- data.frame(rep(NA, 13))
MAG$expr <- mmul_datExpr.norm[which(rownames(mmul_datExpr.norm) == rownames(mmul_HC.fullModel)[which(mmul_HC.fullModel$external_gene_name == 'MAG')]),]
MAG[,1] <- MAG$expr/mean(MAG$expr[which(HC_datMeta$MIA == 'CON')])
colnames(MAG) <- c('norm', 'expr')

# combined the expression of all genes in to one dataframe
myelin_expr <- data.frame(HC_datMeta$MIA)
colnames(myelin_expr) <- c('Group')
myelin_expr$OLIG2 <- OLIG2$norm
myelin_expr$SOX10 <- SOX10$norm
myelin_expr$MYRF <- MYRF$norm
myelin_expr$PLLP <- PLLP$norm
myelin_expr$MBP <- MBP$norm
myelin_expr$CNP <- CNP$norm
myelin_expr$MOG <- MOG$norm
myelin_expr$MAG <- MAG$norm
myelin_expr

library(reshape2)
myelin_expr_melt <- melt(myelin_expr)
colnames(myelin_expr_melt) <- c('Group', 'Gene', 'norm_expr')

pvals <- c(paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'OLIG2')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(hg_HC.fullModel$external_gene_name == 'SOX10')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'MYRF')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'PLLP')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'MBP')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'CNP')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'MOG')], 3)),
           paste0('p =\n', round(mmul_HC.fullModel$P.Value[which(mmul_HC.fullModel$external_gene_name == 'MAG')], 3)))  

library(ggpubr)

y.position <- apply(myelin_expr[,c(2:9)], 2, function(x) max(x, na.rm = TRUE)) + 0.05

myelin_expr_plot <- ggplot(myelin_expr_melt, aes(x = Gene, y = norm_expr, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data=data.frame(), aes(x=myelin_genes, y=y.position, label=pvals), size = 2, inherit.aes = FALSE) +
  ggtitle('DE of Select Genes in HC') +
  theme_bw() +
  ylim(0.9, 1.6) + 
  ylab(expression('Centered Expr (log '[2]*' (CPM))')) +
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

plot(myelin_expr_plot)

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     cell_type_heatmaps, specificity_scale_bar, RFC2, DGCR2, MLC1, KDSR,
     PACS1, PGBD3, GRM3, GRM8, TRIO, GRM1, GRIK2, GRID2, myelin_expr_plot,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (8) logFC correlation for all genes across 4 regions ###############################################################

DLPFC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_V1

# order both topTables by row name so all genes are in the same corresponding rows
DLPFC_ordered <- DLPFC_datExpr[order(rownames(DLPFC_datExpr)),]
ACC_ordered <- ACC_datExpr[order(rownames(ACC_datExpr)),]
HC_ordered <- HC_datExpr[order(rownames(HC_datExpr)),]
V1_ordered <- V1_datExpr[order(rownames(V1_datExpr)),]

library(GGally)

pairsDat <- cbind(DLPFC_ordered$logFC,
                  ACC_ordered$logFC,
                  HC_ordered$logFC,
                  V1_ordered$logFC)
rownames(pairsDat) <- rownames(DLPFC_ordered)
colnames(pairsDat) <- c('DLPFC', 'ACC', 'HC', 'V1')
head(pairsDat)

pairsDat <- data.frame(pairsDat)

all_genes_logFC_cor_plot <- ggpairs(pairsDat,
                          upper = list(continuous = wrap("cor", size = 3, color = 'black')),
                          lower = list(continuous = wrap("points", alpha = 1, size = 0.3))) +
  ggtitle("All Genes logFC Correlation") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin = unit(c(10, 5, 5, 5), 'mm'))

all_genes_logFC_cor_plot <- ggmatrix_gtable(all_genes_logFC_cor_plot)

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     cell_type_heatmaps, specificity_scale_bar, RFC2, DGCR2, MLC1, KDSR,
     PACS1, PGBD3, GRM3, GRM8, TRIO, GRM1, GRIK2, GRID2, myelin_expr_plot, 
     all_genes_logFC_cor_plot,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')



### (9) logFC correlation for genes with p < 0.005 across 4 regions ####################################################

DLPFC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
ACC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_ACC
HC_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_HC
V1_datExpr <- DGE_result_topTables$mmul_CON_vs_MIA_V1


# order both topTables by row name so all genes are in the same corresponding rows
DLPFC_ordered <- DLPFC_datExpr[order(rownames(DLPFC_datExpr)),]
ACC_ordered <- ACC_datExpr[order(rownames(ACC_datExpr)),]
HC_ordered <- HC_datExpr[order(rownames(HC_datExpr)),]
V1_ordered <- V1_datExpr[order(rownames(V1_datExpr)),]

genelist <- rownames(DLPFC_ordered)[which(DLPFC_ordered$P.Value < 0.005)]
genelist <- append(genelist, rownames(ACC_ordered)[which(ACC_ordered$P.Value < 0.005)])
genelist <- append(genelist, rownames(HC_ordered)[which(HC_ordered$P.Value < 0.005)])
genelist <- append(genelist, rownames(V1_ordered)[which(V1_ordered$P.Value < 0.005)])

genelist <- unique(genelist)

# order both topTables by row name so all genes are in the same corresponding rows
pairsDat <- cbind(DLPFC_ordered$logFC[which(rownames(DLPFC_ordered) %in% genelist)],
                  ACC_ordered$logFC[which(rownames(ACC_ordered) %in% genelist)],
                  HC_ordered$logFC[which(rownames(HC_ordered) %in% genelist)],
                  V1_ordered$logFC[which(rownames(V1_ordered) %in% genelist)])
rownames(pairsDat) <- rownames(DLPFC_ordered)[which(rownames(DLPFC_ordered) %in% genelist)]
colnames(pairsDat) <- c('DLPFC', 'ACC', 'HC', 'V1')
head(pairsDat)

pairsDat <- data.frame(pairsDat)

DE_genes_logFC_cor_plot <- ggpairs(pairsDat,
                          upper = list(continuous = wrap("cor", size = 3, color = 'black')),
                          lower = list(continuous = wrap("points", alpha = 1, size = 0.3))) +
  ggtitle("p < 0.005 Genes logFC Correlation") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin = unit(c(2, 5, 2, 5), 'mm'))

DE_genes_logFC_cor_plot <- ggmatrix_gtable(DE_genes_logFC_cor_plot)

save(regional_volcano_plot, num_DE_genes_plot, regional_GO_enrichment_plot,
     cell_type_heatmaps, specificity_scale_bar, RFC2, DGCR2, MLC1, KDSR,
     PACS1, PGBD3, GRM3, GRM8, TRIO, GRM1, GRIK2, GRID2, myelin_expr_plot,
     all_genes_logFC_cor_plot, DE_genes_logFC_cor_plot, regional_EWCE_enrichment_plot,
     file = './working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')


