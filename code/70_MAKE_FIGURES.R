########################################################################################################################
#
#   70_MAKE_FIGURES.R
#
#   Make figures from graphs saved during the analysis and export them as PDFs
# 
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
########################################################################################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

library(gridExtra)
library(ggplot2)
library(Cairo)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')



### (1) Figure 1 #######################################################################################################

source('./code/11_QC_AND_EDA_FIGURES.R')
load('./working_data/plots/11_QC_AND_EDA_FIGURES.RData')

g1 <- ggplot()
plot_list <- list(g1, g1, g1, ePC_plot)

layout_matrix <- rbind(c(1, 1, 1, 1, 2, 2, 2),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(3, 3, 3, 3, 2, 2, 2),
                       c(3, 3, 3, 3, 4, 4, 4),
                       c(3, 3, 3, 3, 4, 4, 4))

cairo_pdf('./figures/Fig 1 Draft.pdf', width = 5.5, height = 4)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (2) Figure S1 ######################################################################################################

load('./working_data/plots/11_QC_AND_EDA_FIGURES.RData')

plot_list <- list(trait_cor_plot, pre_norm_density, post_norm_density, outlier_plot)

layout_matrix <- rbind(c(1, 1, 1),
                       c(1, 1, 1),
                       c(1, 1, 1),
                       c(2, 3, 4))


cairo_pdf('./figures/Fig S1 Draft.pdf', width = 6, height = 8)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (3) Figure S2 ######################################################################################################

source('./code/21_LIMMA-VOOM_DGE_FIGURES.R')
load('./working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(combined_volcano_plot, combined_GO_enrichment_plot, combined_EWCE_enrichment_plot)

layout_matrix <- rbind(c(1, 1, 2),
                       c(1, 1, 3))


cairo_pdf('./figures/Fig S2 Draft.pdf', width = 7.5, height = 5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (4) Figure S3 ######################################################################################################

load('./working_data/plots/21_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(poly_volcano_plot, poly_logFC_plot, poly_GO_enrichment_plot)

layout_matrix <- rbind(c(1, 1, 1, 1),
                       c(1, 1, 1, 1),
                       c(2, 2, 3, 3),
                       c(2, 2, 3, 3))


cairo_pdf('./figures/Fig S3 Draft.pdf', width = 7.5, height = 7.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (5) Figure 2 #######################################################################################################

source('./code/22_LIMMA-VOOM_DGE_FIGURES.R')
load('./working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(regional_volcano_plot, all_genes_logFC_cor_plot, g1, g1, num_DE_genes_plot)

layout_matrix <- rbind(c(1, 1, 1, 1, 5, 5, 5),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(3, 3, 4, 4, 4, 4, 4),
                       c(3, 3, 4, 4, 4, 4, 4))


cairo_pdf('./figures/Fig 2 Draft Pt 1.pdf', width = 7.5, height = 6.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()

###

source('./code/31_WGCNA_FIGURES.R')
load('./working_data/plots/31_WGCNA_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, g1, g1, ME_trait_associations, g1)

layout_matrix <- rbind(c(1, 1, 1, 1, 5, 5, 5),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(1, 1, 1, 1, 2, 2, 2),
                       c(3, 3, 4, 4, 4, 4, 4),
                       c(3, 3, 4, 4, 4, 4, 4))


cairo_pdf('./figures/Fig 2 Draft Pt 3.pdf', width = 7.5, height = 6.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (5) Figure S4 ######################################################################################################

source('./code/22_LIMMA-VOOM_DGE_FIGURES.R')
load('./working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(DE_genes_logFC_cor_plot, g1, g1, g1)

layout_matrix <- rbind(c(1, 2),
                       c(3, 4))


cairo_pdf('./figures/Fig S4 Draft Pt 1.pdf', width = 7.5, height = 5.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()

###

source('./code/31_WGCNA_FIGURES.R')
load('./working_data/plots/31_WGCNA_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, ePC_reg_plot, g1, g1)

layout_matrix <- rbind(c(1, 2),
                       c(3, 4))

cairo_pdf('./figures/Fig S4 Draft Pt 2.pdf', width = 7.5, height = 5.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (5) Figure S6 ######################################################################################################

source('./code/22_LIMMA-VOOM_DGE_FIGURES.R')
load('./working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(regional_EWCE_enrichment_plot, regional_GO_enrichment_plot, g1, g1)

layout_matrix <- rbind(c(3, 2),
                       c(1, 2),
                       c(1, 2),
                       c(1, 2),
                       c(1, 2),
                       c(4, 2))


cairo_pdf('./figures/Fig S6 Draft.pdf', width = 7.5, height = 4.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (6) Figure 3 Pt 1 ##################################################################################################

source('./code/31_WGCNA_FIGURES.R')
load('./working_data/plots/31_WGCNA_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, MEsteelblue_boxplot, MEsteelblue_hubs, MEsteelblue_GO)

layout_matrix <- rbind(c(1, 2),
                       c(3, 4))


cairo_pdf('./figures/Fig 3 Draft Pt 1.pdf', width = 3.5, height = 3.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (6) Figure 3 Pt 2 ##################################################################################################

source('./code/32_WGCNA_FIGURES.R')
load('./working_data/plots/32_WGCNA_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(combined_ME_barplot, g1, g1)

layout_matrix <- rbind(c(1, 1),
                       c(3, 4))


cairo_pdf('./figures/Fig 3 Draft Pt 2.pdf', width = 3.5, height = 3.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (7) Figure 5 Pt 1 ##################################################################################################

source('./code/22_LIMMA-VOOM_DGE_FIGURES.R')
load('./working_data/plots/22_LIMMA-VOOM_DGE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, myelin_expr_plot, g1, g1, g1, g1, g1, g1)

layout_matrix <- rbind(c(1, 1, 2, 2, 2, 2, 2, 2),
                       c(3, 3, 4, 4, 4, 5, 5, 5),
                       c(6, 6, 7, 7, 7, 8, 8, 8))


cairo_pdf('./figures/Fig 5 Draft Pt 1.pdf', width = 5.5, height = 5.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (8) Figure 5 Pt 2 ##################################################################################################

source('./code/32_WGCNA_FIGURES.R')
load('./working_data/plots/32_WGCNA_FIGURES.RData')

g1 <- ggplot()

module_GO_plots <- plot_grid(MEturquoise_GO, MEpaleturquoise_GO, nrow = 2, align = 'hv')

plot_list <- list(g1, g1, 
                  MEturquoise_boxplot, MEturquoise_hubs, module_GO_plots, 
                  MEpaleturquoise_boxplot, MEpaleturquoise_hubs)

layout_matrix <- rbind(c(1, 1, 2, 2, 2, 2, 2, 2),
                       c(3, 3, 4, 4, 4, 5, 5, 5),
                       c(6, 6, 7, 7, 7, 5, 5, 5))


cairo_pdf('./figures/Fig 5 Draft Pt 2.pdf', width = 5.5, height = 5.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (9) Figure S7 ######################################################################################################

source('./code/33_WGCNA_FIGURES.R')
load('./working_data/plots/33_WGCNA_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(all_module_cell_type_enrichment, all_modules_GO)

layout_matrix <- rbind(c(1, 1, 2, 2, 2))


cairo_pdf('./figures/Fig S7 Draft.pdf', width = 7.5, height = 9.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (10) Figure 4 ######################################################################################################

source('./code/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.R')
load('./working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')

g1 <- ggplot()


plot_list <- list(PIWIL2, g1, g1, TE_regional_volcano_plot)

layout_matrix <- rbind(c(1, 2, 2, 3),
                       c(1, 2, 2, 3),
                       c(4, 4, 4, 4),
                       c(4, 4, 4, 4),
                       c(4, 4, 4, 4))


cairo_pdf('./figures/Fig 4 Draft.pdf', width = 7.5, height = 4.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (11) Figure S8 #####################################################################################################

source('./code/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.R')
load('./working_data/plots/41_TRANSPOSABLE_ELEMENTS_DE_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(TE_ePC_plot, TE_pre_norm_density, TE_post_norm_density, 
                  TE_combined_volcano_plot, PGBD3, PGBD5)

layout_matrix <- rbind(c(1, 1, 1, 2, 2, 3, 3),
                       c(4, 4, 4, 5, 5, 6, 6))


cairo_pdf('./figures/Fig S8 Draft.pdf', width = 7, height = 5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (12) Figure 6 ######################################################################################################

source('./code/50_DISEASE_OVERLAP_FIGURES.R')
load('./working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(all_RRHO,
                  DE_disease_enrichment)

layout_matrix <- rbind(c(1),
                       c(2))


cairo_pdf('./figures/Fig 6 Draft Pt 1.pdf', width = 3.5, height = 3.5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (13) Figure S9 #####################################################################################################

source('./code/50_DISEASE_OVERLAP_FIGURES.R')
load('./working_data/plots/50_DISEASE_OVERLAP_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, DLPFC_RRHO, g1,
                  module_disease_enrichment)

layout_matrix <- rbind(c(1, 2, 2, 2, 2, 3),
                       c(4, 4, 4, 4, 4, 4),
                       c(4, 4, 4, 4, 4, 4))


cairo_pdf('./figures/Fig S9 Draft.pdf', width = 5.5, height = 6)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (13) Figure S10 ####################################################################################################

source('./code/60_STEREOTYPIES_FIGURES.R')
load('./working_data/plots/60_STEREOTYPIES_FIGURES.RData')

g1 <- ggplot()

plot_list <- list(g1, stereo_plot, MEdarkorange_stereo_plot, g1, stereo_cor_mat)

layout_matrix <- rbind(c(1, 2, 2, 2, 1, 1, 3, 3, 3, 4),
                       c(1, 2, 2, 2, 1, 1, 3, 3, 3, 4),
                       c(1, 2, 2, 2, 1, 1, 3, 3, 3, 4),
                       c(1, 2, 2, 2, 1, 1, 3, 3, 3, 4),
                       c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5))


cairo_pdf('./figures/Fig S10 Draft.pdf', width = 6.5, height = 5)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()



### (14) Combine final figures #########################################################################################

library(pdftools)

setwd('C:/Users/npage/Github/CONTE_NHP_MIA')

figure_list <- c('./figures/final/Fig 1.pdf',
                 './figures/final/Fig 2.pdf',
                 './figures/final/Fig 3.pdf',
                 './figures/final/Fig 4.pdf',
                 './figures/final/Fig 5.pdf',
                 './figures/final/Fig 6.pdf',
                 './figures/final/Fig S1.pdf',
                 './figures/final/Fig S2.pdf',
                 './figures/final/Fig S3.pdf',
                 './figures/final/Fig S4.pdf',
                 './figures/final/Fig S5.pdf',
                 './figures/final/Fig S6.pdf',
                 './figures/final/Fig S7.pdf',
                 './figures/final/Fig S8.pdf',
                 './figures/final/Fig S9.pdf',
                 './figures/final/Fig S10.pdf')

pdf_combine(figure_list, output = './figures/final/Conte NHP MIA Figures.pdf')












