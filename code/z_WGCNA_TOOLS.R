#####################################################################################################################################
#
#   z_WGCNA_TOOLS.R 
#
#   Contains home-made functions that are useful for running and streamlining WGCNA.
#
#   Nicholas Page, July 2019
#   Gandal and Geschwind Labs, UCLA
#
#####################################################################################################################################

library(gProfileR)
library(ggplot2)
library(WGCNA)
library(pSI)

# Regresses the covariates specified by the user out of the expression matrix, 
# it is called as part of the function runCovariateRegression
regressCovariates <- function(datExpr, datMeta, model, to_regress) {

  beta <- (solve(t(model)%*%model)%*%t(model))%*%t(datExpr)
  
  # regress out beta values from data
  to_regress <- (as.matrix(model[,which(colnames(model) %in% to_regress)]) %*% (as.matrix(beta[which(colnames(model) %in% to_regress),])))
  datExpr.reg <- datExpr - t(to_regress)
  
  ## This datExpr.reg is now a technical variable corrected matrix.
  rownames(datExpr.reg) <- rownames(datExpr)
  colnames(datExpr.reg) <- colnames(datExpr)
  
  return(datExpr.reg)
}

# Correlation panel function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if ((class(x) == "numeric" | class(x) == "integer") & (class(y) == "numeric" | class(y) == "integer")) {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.2/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Check relationship between top expression PCs and Covariates before and after Regression
runCovariateRegression <- function(datExpr, datMeta, model, all_covariates, to_regress, figureDir, figureName) {
  
  ePCs <- prcomp(t(scale(t(datExpr), scale = TRUE)), center = FALSE)
  ePCs <- data.frame(ePCs$rotation)
  
  pdf(paste0(figureDir, figureName, '_unregressed_trait_correlations.pdf'), 24, 20)
  pairs(cbind(all_covariates, ePCs[,c(1:5)]), pch = 19, upper.panel = panel.cor, main = paste("Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values","\n",  sep = ""),cex.labels=1.2)
  dev.off()
  
  ## Regress out covariates
  datExpr.reg <- regressCovariates(datExpr, datMeta, model, to_regress)
  
  ePCs <- prcomp(t(scale(t(datExpr.reg), scale = TRUE)), center = FALSE)
  ePCs <- data.frame(ePCs$rotation)
  
  pdf(paste0(figureDir, figureName, '_regressed_trait_correlations.pdf'), 24, 20)
  pairs(cbind(all_covariates, ePCs[,c(1:5)]), pch = 19, upper.panel = panel.cor, main = paste("Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values","\n",  sep = ""),cex.labels=1.2)
  dev.off()
  
  return(datExpr.reg)
}

# Relate Dendrogram to traits (this script needs to be adapted for each individual dataset/covariates...)
relateDendroTraits <- function(geneTree, TOM, datMeta, datExpr, outputName) {
  
  geneSigs=matrix(NA,nrow = 3,ncol=nrow(datExpr)) # create a vector to hold the data
  
  for( i in 1:ncol(geneSigs)) {
    # calculate r correlation value for numeric variables
    # calculate adjusted R^2s square-root for categorical variables (factors)
    exprvec=as.numeric(datExpr[i,]) # get the expression vector for ith gene
    
    MIA_r=sqrt(max(summary(lm(exprvec~as.factor(datMeta$MIA=="MIA")))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.factor(datMeta$MIA=="MIA")))$coefficients[2,1])
    poly1_r=sqrt(max(summary(lm(exprvec~as.factor(datMeta$Group=="poly1")))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.factor(datMeta$Group=="poly1")))$coefficients[2,1])
    poly2_r=sqrt(max(summary(lm(exprvec~as.factor(datMeta$Group=="poly2")))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.factor(datMeta$Group=="poly2")))$coefficients[2,1])
    #DLPFC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="DLPFC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="DLPFC"))))$coefficients[2,1]) # for DLPFC
    #ACC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="ACC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="ACC"))))$coefficients[2,1]) # for ACC
    #HC_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="HC"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="HC"))))$coefficients[2,1]) # for HC
    #V1_r <- sqrt(max(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="V1"))))$adj.r.squared,0)) * sign(summary(lm(exprvec~as.numeric(as.factor(datMeta$Region=="V1"))))$coefficients[2,1]) # for V1
    geneSigs[,i]=c(MIA_r,poly1_r,poly2_r)
  }
    
  colnames(geneSigs)=rownames(datExpr)
  rownames(geneSigs)=c("MIA", "poly1", "poly2")
  
  geneSigsColor=matrix(NA,nrow=nrow(geneSigs),ncol=nrow(datExpr)) # create a vector to hold the data
  for ( i in 1:nrow(geneSigsColor)) {
    geneSigsColor[i,] = numbers2colors(as.numeric(geneSigs[i,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 
  }
  
  rownames(geneSigsColor)=rownames(geneSigs)
  colnames(geneSigsColor)=colnames(geneSigs)
  
  # Try out tree cutting parameters
  mColorh <- mLabelh <- colorLabels <- NULL  
  for (minModSize in c(50,100,150)) {
    for (dthresh in c(0.1,0.2,0.25)) {
      for (ds in c(2,4)) {
        print("Trying parameters:")
        print(c(minModSize,dthresh,ds))
        tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE, minClusterSize = minModSize, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(1-TOM))
        merged <- mergeCloseModules(exprData = t(datExpr),colors = tree$labels,
                                    cutHeight = dthresh)
        mColorh <- cbind(mColorh,labels2colors(merged$colors))
        mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
      }
    }
  } 
  mColorh1=cbind(mColorh,t(geneSigsColor))
  mLabelh1=c(mLabelh,rownames(geneSigsColor))
  save(mColorh1,mLabelh1,geneSigs,geneSigsColor,file=paste0('./', outputName,".RData"))
  
  pdf(paste0('./', outputName,".pdf"),height=25,width=20)
  plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
  dev.off()
  
  return(geneSigsColor)
} 



finalDendroCut <- function(datExpr,geneTree,dissTOM,geneSigsColor,softPower,processedDataDir,mms,ds,dthresh) {
  if (!file.exists(paste0(processedDataDir,"modules.RData"))) {
    tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
    
    merged <- mergeCloseModules(exprData = t(datExpr), colors = tree$labels, cutHeight = dthresh)
    mColorh <- cbind(labels2colors(merged$colors),t(geneSigsColor))
    mLabelh <- c("Merged Colors",rownames(geneSigsColor))
    
    pdf(paste0(processedDataDir,"_Signed_Dendro_Parameters.pdf"),height=10,width=16)
    plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh));
    dev.off()
    
    mergedColors = labels2colors(merged$colors);
    
    # Eigengenes of the new merged modules:
    MEList=moduleEigengenes(t(datExpr), colors = mergedColors, softPower = softPower, nPC=1)
    MEs=MEList$eigengenes
    MEs=orderMEs(MEs)
    rownames(MEs) = colnames(datExpr)
    moduleColors = mergedColors
    # names(moduleColors) <- rownames(datExpr)
    save(geneTree,moduleColors,MEs,file=paste0(processedDataDir,"_Modules.RData"))
  }
}



findEnrichment <- function(MEs, datMeta, datProbes, moduleColors) {
  
  #Load cell-type expression signatures for enrichment
  zhang.datExpr = read.csv("./datExpr.zhangHuman.avgForPSI.csv",row.names=1)[,-1]
  pSI.output = specificity.index(pSI.in=zhang.datExpr,bts=100,p_max=.1, e_min=0.3)
  pSI.count(pSI.output)
  
  # fits a linear model to the data that will be used to determine module enrichment
  fit <- lmFit(t(MEs), mod)
  fit_contrasts <- contrasts.fit(fit, contrast.matrix)
  results <- eBayes(fit_contrasts)
  
  # get the results of the linear model
  tt <- topTable(results, number = Inf, coef = 1)
  DLPFC_tt <- topTable(results, number = Inf, coef = 2)
  ACC_tt <- topTable(results, number = Inf, coef = 3)
  HC_tt <- topTable(results, number = Inf, coef = 4)
  V1_tt <- topTable(results, number = Inf, coef = 5)
  
  
  pdf(file = './30_WGCNA_Module_Enrichments.pdf', width = 14, height = 12)
  
  for(i in 1:dim(MEs)[2]) {
    
    # get the module eigengene
    me = MEs[,i]
    m = colnames(MEs)[i]
    moduleColor = c = substr(m,3,nchar(m))
    
    moduleGenes = datProbes$ensembl_gene_id[which(moduleColors == moduleColor)]
    
    dat= cbind(data.frame(ME = me), datMeta)
    
    g1=ggplot(dat, aes(x=MIA, y=ME)) + geom_boxplot() + 
      geom_point(aes(color=Group), position=position_jitter(.1),size=2) + 
      facet_wrap(~Region) + ggtitle(m, subtitle = paste0('MIA: p=',tt$P.Value[which(rownames(tt) == m)])) + 
      theme(plot.title  = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5)) + xlab("")
    
    
    pTable <- data.frame("P-Value" = c(tt$P.Value[which(rownames(tt) == m)],
                                       DLPFC_tt$P.Value[which(rownames(DLPFC_tt) == m)],
                                       ACC_tt$P.Value[which(rownames(ACC_tt) == m)],
                                       HC_tt$P.Value[which(rownames(HC_tt) == m)],
                                       V1_tt$P.Value[which(rownames(V1_tt) == m)]))
    rownames(pTable) <- c('ALL','DLPFC', 'ACC', 'HC', 'V1')
    g.pTable = tableGrob(pTable)

    kME = signedKME(t(datExpr.norm.reg), MEs,corFnc = "bicor")
    
    #Gene Ontology
    idx=which(moduleColors == moduleColor)
    query = rownames(kME)[idx[order(kME[idx,paste0('kME',c)],decreasing = T)]]
    go = gprofiler(query, organism="mmulatta", custom_bg = datProbes$ensembl_gene_id, 
                   correction_method = "fdr",hier_filtering = "strong", ordered_query = T, significant = T, exclude_iea = F,
                   region_query = F,max_p_value = 0.05, max_set_size=1000, numeric_ns = "",
                   include_graph = F,src_filter = c("GO", "KEGG"))
    go = go[order(go$p.value)[1:min(10,nrow(go))],]
    g3=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")+theme(axis.text.y = element_text(size=8)) + ggtitle("Gene Ontolgoy")
    
    
    #Hub Gene Network
    hubGenes = moduleGenes[order(kME[moduleGenes,i], decreasing = T)[1:20]]
    hubGene.symbols = datProbes$external_gene_name[match(hubGenes, datProbes$ensembl_gene_id)]
    adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "bicor", power=10)
    adjMat[adjMat < quantile(adjMat,0.1)]=0
    graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
    plotcord= data.frame(layout_with_fr(graph))
    colnames(plotcord) = c("X1","X2")
    edgelist <- get.edgelist(graph,names = F)
    edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
    colnames(edges) <- c("X1","Y1","X2","Y2")
    plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
    g4=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
      geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),data=plotcord) +
      theme_classic() + labs(x="", y="") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    
    moduleGenes.hs = datProbes$hsapiens_gene_id[match(moduleGenes, datProbes$ensembl_gene_id)]
    moduleGenes.hs <- moduleGenes.hs[which(moduleGenes.hs != '')]
    f = fisher.iteration(pSI.output, moduleGenes.hs,p.adjust = T)
    f$CellType = rownames(f)
    f$log10fdr = -log10(f[,1])
    f$Dataset = "Zhang"
    f = f[,-c(1:4)]
    
    g5=ggplot(f,aes(x=CellType,y=log10fdr, fill=Dataset)) + geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept = -log10(0.05),lty=2) + xlab("") + ggtitle("Cell-Type Enrichment") +theme(axis.text.y = element_text(size=6))
    
    grid.arrange(grobs=list(g1,g3,g4,g5,g.pTable), layout_matrix=rbind(c(1,1, 4, 4,2,2 ),c(1,1, 3,3,3,3),c(5,5,3,3, 3,3)),padding=unit(2, "line"))
    
  }
  
  dev.off()
}










