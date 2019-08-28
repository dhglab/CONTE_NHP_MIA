########################################################################################################################
#
#   z_PREPARE_SPECIFICITY_FOR_EWCE.R 
#
#   Load scRNA-seq data compiled by the PsychEncode Consortium (PEC) and generate a cell-type specificity for the 
#   expression of each gene to determine cell-type enrichment in gene sets.
#
#   Nicholas Page, August 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Setup Lake + PsychENCODE (2018) data
#   (2) Calculate specificity metric for Lake + PsychENCODE (2018) data
#
########################################################################################################################

### (1) Setup Lake + PsychENCODE (2018) data ###########################################################################

PEC_sc_expr = read.table('./raw_data/DER-22_Single_cell_expression_raw_UMI.tsv', sep = '\t') # 17176 genes, 27412 cells
PEC_sc_annot = data.frame('cell_id' = colnames(PEC_sc_expr))

# get broad and specific cell types from the scRNA-seq data
PEC_sc_annot$specific_celltype = gsub('\\..*', '', PEC_sc_annot$cell_id)
PEC_sc_annot$broad_celltype = PEC_sc_annot$specific_celltype
PEC_sc_annot$broad_celltype[grep('Ex', PEC_sc_annot$broad_celltype)] = 'ExN'
PEC_sc_annot$broad_celltype[grep('In', PEC_sc_annot$broad_celltype)] = 'InN'
PEC_sc_annot$broad_celltype[grep('Microglia', PEC_sc_annot$broad_celltype)] = 'Micro'

# remove cells with an unkown cell-type
idx = which(PEC_sc_annot$specific_celltype == 'NA')
PEC_sc_annot = PEC_sc_annot[-idx,]
PEC_sc_expr = PEC_sc_expr[,-idx] # 17176 genes, 27380 cells

save(PEC_sc_expr, PEC_sc_annot, file = './working_data/PEC_Lake_AdultNeural-scRNAseq_for_EWCE.RData')



### (2) Calculate specificity metric for Lake + PsychENCODE (2018) data ################################################

load('./working_data/PEC_Lake_AdultNeural-scRNAseq_for_EWCE.RData')

PEC_sc_expr = as.matrix(PEC_sc_expr) # 17176 genes, 27380 cells

# remove uninformative genes
PEC_sc_expr_DROPPED = drop.uninformative.genes(exp = PEC_sc_expr, level2annot = PEC_sc_annot$broad_celltype) # 17174 genes, 27380 cells
annotLevels = list(specific_celltype = PEC_sc_annot$specific_celltype, broad_celltype = PEC_sc_annot$broad_celltype)

# generate cell-type specificity matrix
fileNames_PEC_Lake_Adult_Human_nucSeq = generate.celltype.data(exp = PEC_sc_expr_DROPPED, 
                                                               annotLevels = annotLevels, 
                                                               groupName = "PEC_Lake_Adult_Human_nucSeq")
print(fileNames_PEC_Lake_Adult_Human_nucSeq)

# Proceed to EWCE
