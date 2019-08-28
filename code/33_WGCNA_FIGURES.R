########################################################################################################################
#
#   33_WGCNA_FIGURES.R 
#
#   Generate figures for WGCNA on 4 y/o non-human primate RNA-seq 
#   samples from DLPFC, ACC, HC, and V1
#
#   Nicholas Page, August 2019
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
#   (1) Cell-type enrichment of all WGCNA modules
#   (2) Top GO term enrichment for all WGCNA modules
#
########################################################################################################################

### (1) Cell-type enrichment of all WGCNA modules ######################################################################

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

background_adult = rownames(ctd[[level1]]$specificity)
length(background_adult)
background_adult = background_adult[which(background_adult %in% mmul_CON_vs_MIA_ALL$external_gene_name)]
length(background_adult)

ewce_results <- data.frame()

for(i in 1:dim(MEs)[2]) {
  
  ME_i = names(moduleColors[which(moduleColors == str_replace(colnames(MEs)[i], 'ME', ''))])
  ME_i = ME_i[which(ME_i %in% background_adult)]
  ME_i_results = bootstrap.enrichment.test(sct_data=ctd, hits=ME_i, bg=background_adult,reps=reps,annotLevel=level1,genelistSpecies="human",sctSpecies="human")
  
  ewce_results <- rbind(ewce_results,
                        data.frame(ME_i_results$results[,c(1,3,5)], 'Module' = colnames(MEs)[i]))
  
}

ewce_results$FDR <- p.adjust(ewce_results$p, method = 'fdr')
ewce_results$sd_from_mean <- ifelse(ewce_results$sd_from_mean > 0, ewce_results$sd_from_mean, 0)
print(ewce_results)

ewce_results$sig <- ifelse(ewce_results$FDR < 0.05, '*', '')

ewce_results$Module <- factor(str_replace(as.character(ewce_results$Module), 'ME', ''),
                             levels = str_replace(as.character(ewce_results$Module), 'ME', ''),
                             labels = str_replace(as.character(ewce_results$Module), 'ME', ''))

all_module_cell_type_enrichment <- ggplot(data = ewce_results, aes(x = ewce_results$CellType, y = ewce_results$sd_from_mean, 
                                                                   fill = ewce_results$Module, 
                                                                   color = ewce_results$Module)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(facets = vars(Module), nrow = 10, ncol = 3) +
  geom_text(aes(label=ewce_results$sig), vjust=0.5, color="red",
            position = position_dodge(0), size = 8) +
  scale_color_manual(values = unique(str_replace(as.character(ewce_results$Module), 'ME', ''))) +
  scale_fill_manual(values = unique(str_replace(as.character(ewce_results$Module), 'ME', ''))) +
  ggtitle('Cell-type Enrichment from PsychEncode') +
  ylab('Stdev from Mean') +
  ylim(c(0, 35)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 1, margin=margin(-5,0,0,0)),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text =  element_text(size = 6),
        plot.margin=unit(c(2,2,2,2),"mm"))

plot(all_module_cell_type_enrichment)

save(all_module_cell_type_enrichment,
     file = './working_data/plots/33_WGCNA_FIGURES.RData')



### (2) Top GO term enrichment for all WGCNA modules ###################################################################

source('./code/z_WGCNA_TOOLS.R')

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

datExpr.norm.reg <- runCovariateRegression(datExpr = datExpr.norm, datMeta = datMeta, model = mod,
                                           all_covariates = datMeta[,c(2:4,6:12)], to_regress = c('MG_seqPC1', 'MG_seqPC2', 'MG_seqPC3', 'MG_seqPC4'), 
                                           figureDir = './', figureName = 'mmul_gene_WGCNA')

go_results <- data.frame()

for(j in 1:dim(MEs)[2]) {
  
  # get the module eigengene
  me = MEs[,j]
  m = colnames(MEs)[j]
  moduleColor = c = substr(m,3,nchar(m))
  
  moduleGenes = datProbes$ensembl_gene_id[which(moduleColors == moduleColor)]
  
  dat= cbind(data.frame(ME = me), datMeta)
  
  kME = signedKME(t(datExpr.norm.reg), MEs,corFnc = "bicor")
  
  #Gene Ontology
  idx=which(moduleColors == moduleColor)
  query = rownames(kME)[idx[order(kME[idx,paste0('kME',c)],decreasing = T)]]
  go = gprofiler(query, organism="mmulatta", custom_bg = datProbes$ensembl_gene_id, 
                 correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = F, exclude_iea = F,
                 region_query = F,max_p_value = 1, max_set_size=1000, numeric_ns = "",
                 include_graph = F,src_filter = c("GO", "KEGG"))
  go = go[order(go$p.value)[1:min(2,nrow(go))],]
  
  go_results <- rbind(go_results,
                      data.frame(go, 'Module' = colnames(MEs)[j]))
  
}

go_results$p.value <- as.numeric(go_results$p.value)
GO_combined <- data.frame(go_results)
GO_combined$term.name <- factor(GO_combined$term.name,
                                   levels = GO_combined$term.name,
                                   labels = GO_combined$term.name)
GO_combined$Module <- factor(str_replace(as.character(GO_combined$Module), 'ME', ''),
                                levels = str_replace(as.character(GO_combined$Module), 'ME', ''),
                                labels = str_replace(as.character(GO_combined$Module), 'ME', ''))

all_modules_GO <- ggplot(data = GO_combined, aes(x = GO_combined$term.name, y = -log10(as.numeric(GO_combined$p.value)), 
                                        color = GO_combined$Module,
                                        fill = GO_combined$Module)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_x_discrete(limits = rev(GO_combined$term.name)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30), limits = c(0, 30)) +
  scale_color_manual(values = unique(str_replace(as.character(GO_combined$Module), 'ME', ''))) +
  scale_fill_manual(values = unique(str_replace(as.character(GO_combined$Module), 'ME', ''))) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +
  ggtitle('GO Enrichment') +
  ylab(expression('\n-log'[10]*'(P-Value)')) +
  theme_bw() +
  #guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 5, l = 0)),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, lineheight = 0.6),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = 'bottom',
        strip.text =  element_text(size = 10),
        legend.margin = margin(t = 0, r = 60, b = 0, l = 0, unit = 'mm'))

# plot(all_modules_GO)

save(all_module_cell_type_enrichment, all_modules_GO,
     file = './working_data/plots/33_WGCNA_FIGURES.RData')






























