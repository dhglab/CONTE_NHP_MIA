########################################################################################################################
#
#   50_CROSS_DISORDER_OVERLAP.R 
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

########################################################################################################################
#
#   Pipeline:
#
#   (1) Perform gene set enrichment with up and down genelists
#   (2) Perform gene set enrichment with DE modules
#   (3) Plot TE expression density post-normalization
#   (4) RRHO on DLPFC DE genes vs disease
#
########################################################################################################################

### (1) Perform gene set enrichment with up and down genelists #########################################################

genelists <- read.csv('./raw_data/Genelists.csv', header = TRUE, sep = ',')
genelists <- genelists[,c(3:9,11,15)]

cross_disorder_DE <- read.csv(file = './raw_data/aat8127_Table_S1.csv', header = TRUE, sep = ',')

ASD_UP <- cross_disorder_DE[which(cross_disorder_DE$ASD.fdr < 0.01 & cross_disorder_DE$ASD.log2FC > 0),]
genelists$ASD_UP <- c(ASD_UP$gene_name, rep('', (dim(genelists)[1] - length(ASD_UP$gene_name))))

ASD_DOWN <- cross_disorder_DE[which(cross_disorder_DE$ASD.fdr < 0.01 & cross_disorder_DE$ASD.log2FC < 0),]
genelists$ASD_DOWN <- c(ASD_DOWN$gene_name, rep('', (dim(genelists)[1] - length(ASD_DOWN$gene_name))))

SCZ_UP <- cross_disorder_DE[which(cross_disorder_DE$SCZ.fdr < 0.00001 & cross_disorder_DE$SCZ.log2FC > 0),]
genelists$SCZ_UP <- c(SCZ_UP$gene_name, rep('', (dim(genelists)[1] - length(SCZ_UP$gene_name))))

SCZ_DOWN <- cross_disorder_DE[which(cross_disorder_DE$SCZ.fdr < 0.00001 & cross_disorder_DE$SCZ.log2FC < 0),]
genelists$SCZ_DOWN <- c(SCZ_DOWN$gene_name, rep('', (dim(genelists)[1] - length(SCZ_DOWN$gene_name))))

ref_background <- cross_disorder_DE$gene_name

genelists <- genelists[,c(9:13)]

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

source('./code/z_OVER_ENRICHMENT_ANALYSIS.R')

query_list <- list(DLPFC_UP, DLPFC_DOWN, ACC_UP, ACC_DOWN, HC_UP, HC_DOWN, V1_UP, V1_DOWN)
query_name_list <- list('DLPFC_UP', 'DLPFC_DOWN', 'ACC_UP', 'ACC_DOWN', 'HC_UP', 'HC_DOWN', 'V1_UP', 'V1_DOWN')

query_background <- mmul_CON_vs_MIA_ALL$external_gene_name

ORA_summary <- list()

mapply(function(query, query_name) {
  
  for(i in 1:dim(genelists)[2]){
    
    term <- genelists[i]
    term <- term[which(term != ''),]
    enrichment <- ORA(query, term, query_background, ref_background)
    ORA_summary <<- rbind(ORA_summary, c('query' = query_name, 
                                         'term' = colnames(genelists)[i], 
                                         'OR' = as.numeric(enrichment[[1]]), 
                                         'p' = as.numeric(enrichment[[2]]), 
                                         'overlap' = as.numeric(enrichment[[9]])))
    
  }
  
}, query_list, query_name_list)

ORA_summary <- data.frame(ORA_summary)

ORA_summary$p <- as.numeric(ORA_summary$p)
ORA_summary$OR <- as.numeric(ORA_summary$OR)
ORA_summary$term <- unlist(ORA_summary$term)
ORA_summary$query <- unlist(ORA_summary$query)

ORA_summary <- ORA_summary[order(ORA_summary$p, decreasing = FALSE),]

ORA_summary$FDR <- p.adjust(ORA_summary$p, method = 'bonferroni')
ORA_summary

ORA_summary$query <- factor(ORA_summary$query,
                            levels = query_name_list,
                            labels = query_name_list)


library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

DE_disease_enrichment <- ggplot(ORA_summary, aes(x = query, y = term, fill = -log10(p))) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(ORA_summary$FDR < 0.05, '*', ifelse(ORA_summary$p < 0.05, '-', ''))), 
            size = 4, color = 'black', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 4)) +
  ggtitle('Disease Enrichment') +
  labs(fill = expression('-log'[10]*'(p)')) +
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
        legend.title = element_text(size = 8),
        plot.margin=unit(c(2,5,5,2),"mm"),
        legend.box.margin=margin(-10,5,-10,-10))

plot(DE_disease_enrichment)

save(DE_disease_enrichment,
     file = './working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')


### (2) Perform gene set enrichment with DE modules ####################################################################

cross_disorder_modules <- read.csv(file = './raw_data/aat8127_Table_S5.csv', header = TRUE, sep = ',')

load('./working_data/WGCNA/WGCNA_sft_14_Final_Modules.RData')
gene2name <- match(names(moduleColors), mmul_annotEns87$ensembl_gene_id)
names(moduleColors) <- mmul_annotEns87$hsapiens_homolog_associated_gene_name[gene2name]

module_ORA_summary <- data.frame()

for(i in 1:length(levels(factor(moduleColors)))){
  
  for(j in 1:length(levels(factor(cross_disorder_modules$Module)))){
    
    query <- names(moduleColors)[which(moduleColors == levels(factor(moduleColors))[i])]
    query <- query[which(query != '')]
    query_name <- levels(factor(moduleColors))[i]
    
    term <- cross_disorder_modules$gene_name[which(cross_disorder_modules$Module == levels(factor(cross_disorder_modules$Module))[j])]
    term <- term[which(term != '')]
    term_name <- levels(factor(cross_disorder_modules$Module))[j]
    
    enrichment <- ORA(query, term, query_background, ref_background)
    module_ORA_summary <- rbind(module_ORA_summary, data.frame('query' = query_name, 
                                                      'term' = term_name, 
                                                      'OR' = as.numeric(enrichment[[1]]), 
                                                      'p' = as.numeric(enrichment[[2]]), 
                                                      'overlap' = as.numeric(enrichment[[9]])))
    
  }
  
}

#module_ORA_summary <- data.frame(module_ORA_summary)

#module_ORA_summary$p <- as.numeric(module_ORA_summary$p)
#module_ORA_summary$OR <- as.numeric(module_ORA_summary$OR)

module_ORA_summary <- module_ORA_summary[order(module_ORA_summary$p, decreasing = FALSE),]
# module_ORA_summary <- module_ORA_summary[order(module_ORA_summary$OR, decreasing = FALSE),]

module_ORA_summary$FDR <- p.adjust(module_ORA_summary$p, method = 'bonferroni')
module_ORA_summary

overlap_modules <- module_ORA_summary[which(module_ORA_summary$FDR < 0.05),]
overlap_modules

library(RColorBrewer)
palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')

module_disease_enrichment <- ggplot(module_ORA_summary, aes(x = query, y = term, fill = -log10(p))) + 
  geom_tile(color = 'white', size = 0.5) +
  geom_text(aes(label = ifelse(module_ORA_summary$FDR < 0.05, '*', '')), 
            size = 4, color = 'black', fontface = 'bold') +
  scale_fill_gradientn(colours = palette(100), limits = c(0, 90)) +
  ggtitle('PsychEncode Module Overlap') +
  labs(fill = expression('-log'[10]*'(p)')) +
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
        legend.title = element_text(size = 8),
        plot.margin=unit(c(2,5,5,2),"mm"),
        legend.box.margin=margin(-10,5,-10,-10))

plot(module_disease_enrichment)

save(DE_disease_enrichment, module_ORA_summary,
     file = './working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')



### (3) RRHO on DE genes vs disease ####################################################################################

library(RRHO)

# prepare monkey list for RRHO
mmul_all <- DGE_result_topTables$mmul_CON_vs_MIA_ALL
mmul_all <- mmul_all[,c(8,1,4)]
mmul_all <- mmul_all[which(mmul_all$external_gene_name != ''),]
mmul_all$sign_p <- -log10(mmul_all$P.Value) * sign(mmul_all$logFC)
mmul_all <- mmul_all[order(mmul_all$sign_p, decreasing = TRUE), c(1,4)]
mmul_all <- mmul_all[which(duplicated(mmul_all$external_gene_name) == FALSE),]

# prepare SCZ gene list for RRHO
SCZ_all <- cross_disorder_DE[,c(8,20,24)]
SCZ_all$sign_p <- -log10(SCZ_all$SCZ.p.value) * sign(SCZ_all$SCZ.log2FC)
SCZ_all <- SCZ_all[order(SCZ_all$sign_p, decreasing = TRUE), c(1,4)]
SCZ_all <- SCZ_all[which(!(is.na(SCZ_all$sign_p))),]
SCZ_all <- SCZ_all[which(duplicated(SCZ_all$gene_name) == FALSE),]

# prepare ASD gene list for RRHO
ASD_all <- cross_disorder_DE[,c(8,14,18)]
ASD_all$sign_p <- -log10(ASD_all$ASD.p.value) * sign(ASD_all$ASD.log2FC)
ASD_all <- ASD_all[order(ASD_all$sign_p, decreasing = TRUE), c(1,4)]
ASD_all <- ASD_all[which(!(is.na(ASD_all$sign_p))),]
ASD_all <- ASD_all[which(duplicated(ASD_all$gene_name) == FALSE),]


all_SCZ_RRHO.result <- RRHO(mmul_all, SCZ_all, BY = TRUE, alternative = 'enrichment', stepsize = 100, log10.ind = TRUE)
all_ASD_RRHO.result <- RRHO(mmul_all, ASD_all, BY = TRUE, alternative = 'enrichment', stepsize = 100, log10.ind = TRUE)


#plot(p)


library(ggplot2)
library(reshape2)

all_SCZ_df = data.frame(melt(all_SCZ_RRHO.result$hypermat), Group='SCZ')
all_ASD_df = data.frame(melt(all_ASD_RRHO.result$hypermat), Group='ASD')

all_df <- rbind(all_SCZ_df, all_ASD_df)
all_df[mapply(is.na, all_df)] <- 0

all_RRHO <- ggplot(all_df, aes(x=Var1,y=Var2,fill=value)) + geom_tile(color=NA) + 
  ggtitle('Disease MIA Overlap') +
  facet_wrap(~Group, nrow=1) + 
  labs(x='DGE rank in MIA', y="DGE rank in Disease", fill = expression('-log'[10]*'(p)')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5),
        strip.text = element_text(size = 10))

plot(all_RRHO)


save(DE_disease_enrichment, module_ORA_summary, all_RRHO,
     file = './working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')



### (4) RRHO on DLPFC DE genes vs disease ##############################################################################

library(RRHO)

# prepare monkey list for RRHO
mmul_DLPFC <- DGE_result_topTables$mmul_CON_vs_MIA_DLPFC
mmul_DLPFC <- mmul_DLPFC[,c(8,1,4)]
mmul_DLPFC <- mmul_DLPFC[which(mmul_DLPFC$external_gene_name != ''),]
mmul_DLPFC$sign_p <- -log10(mmul_DLPFC$P.Value) * sign(mmul_DLPFC$logFC)
mmul_DLPFC <- mmul_DLPFC[order(mmul_DLPFC$sign_p, decreasing = TRUE), c(1,4)]
mmul_DLPFC <- mmul_DLPFC[which(duplicated(mmul_DLPFC$external_gene_name) == FALSE),]

# prepare SCZ gene list for RRHO
SCZ_all <- cross_disorder_DE[,c(8,20,24)]
SCZ_all$sign_p <- -log10(SCZ_all$SCZ.p.value) * sign(SCZ_all$SCZ.log2FC)
SCZ_all <- SCZ_all[order(SCZ_all$sign_p, decreasing = TRUE), c(1,4)]
SCZ_all <- SCZ_all[which(!(is.na(SCZ_all$sign_p))),]
SCZ_all <- SCZ_all[which(duplicated(SCZ_all$gene_name) == FALSE),]

# prepare ASD gene list for RRHO
ASD_all <- cross_disorder_DE[,c(8,14,18)]
ASD_all$sign_p <- -log10(ASD_all$ASD.p.value) * sign(ASD_all$ASD.log2FC)
ASD_all <- ASD_all[order(ASD_all$sign_p, decreasing = TRUE), c(1,4)]
ASD_all <- ASD_all[which(!(is.na(ASD_all$sign_p))),]
ASD_all <- ASD_all[which(duplicated(ASD_all$gene_name) == FALSE),]


DLPFC_SCZ_RRHO.result <- RRHO(mmul_DLPFC, SCZ_all, BY = TRUE, alternative = 'enrichment', stepsize = 100, log10.ind = TRUE)
DLPFC_ASD_RRHO.result <- RRHO(mmul_DLPFC, ASD_all, BY = TRUE, alternative = 'enrichment', stepsize = 100, log10.ind = TRUE)

library(ggplot2)
library(reshape2)

DLPFC_SCZ_df = data.frame(melt(DLPFC_SCZ_RRHO.result$hypermat), Group='SCZ')
DLPFC_ASD_df = data.frame(melt(DLPFC_ASD_RRHO.result$hypermat), Group='ASD')

DLPFC_df <- rbind(DLPFC_SCZ_df, DLPFC_ASD_df)
DLPFC_df[mapply(is.na, DLPFC_df)] <- 0

DLPFC_RRHO <- ggplot(DLPFC_df, aes(x=Var1,y=Var2,fill=value)) + geom_tile(color=NA) + 
  ggtitle('DLPFC Disease MIA Overlap') +
  facet_wrap(~Group, nrow=1) + 
  labs(x='DGE rank in MIA', y="DGE rank in Disease", fill = expression('-log'[10]*'(p)')) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        plot.title = element_text(size = 10, hjust = 0.5),
        strip.text = element_text(size = 10))

plot(DLPFC_RRHO)


save(DE_disease_enrichment, module_ORA_summary, all_RRHO, DLPFC_RRHO,
     file = './working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')

