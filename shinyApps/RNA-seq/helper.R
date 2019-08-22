# helper.R
load("data/data.RData")

th <- theme_bw() + theme(axis.text.x = element_text(size=10), axis.title.x = element_text(size=12), axis.text.y = element_text(size=10), axis.title.y = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5, size=15))

retrieveGeneExpr <- function(normCounts, meta, columns=c("stage", "somite", "somiteNumber"), gene){
  d <- t(as.data.frame(normCounts[normCounts$gene==gene,-1]))
  stopifnot(identical(row.names(d), meta$sample))
  d <- cbind(d, meta[,columns])
  if(ncol(d) <= length(columns)){
    d$cpm <- 0
    d <- d[,c(4,1:3)]
  }else{ colnames(d)[1] <- "cpm" }
  return(d)
}

genesAvailable <- as.character(normCounts$gene)
genesAvailable <- genesAvailable[order(nchar(genesAvailable))]

#### Gene expression
col.stage <- brewer.pal(n=6, "YlOrRd")
boxplotExpr <- function(group_by=4, colour_by=3, gene="Hoxa1"){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 3){
      p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
    }
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 2){
      p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
    }
  }
  return(p)
}

## summary of DE results
## trios
summaryDEtrios <- function(gene="Hoxa1"){
  if(gene != ""){
  ## average
  i.ii <- DEtriosAll[["somiteIvsII"]][DEtriosAll[["somiteIvsII"]]$genes==gene,]
  ii.iii <- DEtriosAll[["somiteIIvsIII"]][DEtriosAll[["somiteIIvsIII"]]$genes==gene,]
  i.iii <- DEtriosAll[["somiteIvsIII"]][DEtriosAll[["somiteIvsIII"]]$genes==gene,]
  
  ## per-stage
  stage <- list()
  for(s in paste0("stage", c(8,18,21,25,27,35))){
    stage[[s]] <- DEtriosStage[[s]][DEtriosStage[[s]]$genes==gene,]
  }
  
  ## plot
  plots <- list()
  rng <- range(c(i.ii$logFC, i.iii$logFC, ii.iii$logFC, stage[["stage8"]][,2:4], stage[["stage18"]][,2:4],
                 stage[["stage21"]][,2:4], stage[["stage25"]][,2:4], stage[["stage27"]][,2:4], stage[["stage35"]][,2:4]))
  
  # avg
  df <- data.frame(var1=factor(c("SI", "SI", "SII"), levels = c("SI", "SII")), 
                   var2=factor(c("SII", "SIII", "SIII"), levels = c("SIII", "SII")), 
                   fc=c(i.ii$logFC, i.iii$logFC, ii.iii$logFC), 
                   fdr=c(i.ii$FDR, i.iii$FDR, ii.iii$FDR))
  df <- df[order(abs(df$fc)),]
  plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) + geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "white"), size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0)) + 
    scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
    theme_minimal() + coord_fixed() + xlab("average") + ylab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(face=2, size=10), axis.title.x = element_text(face=2, size=12))
  
  # stage
  i=2
  for(s in paste0("stage", c(8,18,21,25,27,35))){
    df <- data.frame(var1=factor(c("SI", "SI", "SII"), levels = c("SI", "SII")), 
                     var2=factor(c("SII", "SIII", "SIII"), levels = c("SIII", "SII")), 
                     fc=c(stage[[s]][,2], stage[[s]][,3], stage[[s]][,4]), fdr=stage[[s]]$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) + geom_tile(color = ifelse(abs(df$fc) > log2(1.5), "black", "white"), size = ifelse(abs(df$fc) > log2(1.5), 0.5, 0)) + 
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
      theme_minimal() + coord_fixed() + xlab(s) + ylab("") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text = element_text(face=2, size=10), axis.title.x = element_text(face=2, size=12))
    i <- i+1
  }
  p <- ggarrange(plotlist = plots, ncol = 1, nrow = 7, common.legend = TRUE, legend = "none")
  return(p)
  }
}

## stages
summaryDEstages <- function(gene="Hoxa1"){
  if(gene != ""){
  ## average
  all <- DEstageAll[DEstageAll$genes==gene,]
  
  ## per-somite
  somite <- list()
  for(s in paste0("somite", c("I","II","III"))){
    somite[[s]] <- DEstageSomite[[s]][DEstageSomite[[s]]$genes==gene,]
  }
  
  ## plot
  plots <- list()
  rng <- range(c(all[2:16], somite[["somiteI"]][2:16], somite[["somiteII"]][2:16], somite[["somiteIII"]][2:16]))
  
  # avg
  df <- data.frame(var1=factor(c(rep("8",5), rep("18",4), rep("21",3), rep("25",2), "27"), levels = c("8", "18", "21", "25", "27")),
                   var2=factor(c("18","21","25","27","35", "21","25","27","35", "25","27","35", "27","35", "35"), levels = c("18", "21", "25", "27", "35")),
                   fc=as.numeric(all[,2:16]),
                   fdr=all$FDR)
  df <- df[order(abs(df$fc)),]
  plots[[1]] <- ggplot(df, aes(var1, var2, fill = fc)) + geom_tile(color = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), "black", "white"), size = ifelse(df$fdr < 0.05 & abs(df$fc) > log2(1.5), 0.5, 0)) + 
    scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
    theme_minimal() + coord_fixed() + xlab("average") + ylab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(face=2, size=10), axis.title.x = element_text(face=2, size=12))
  
  # somite
  i=2
  for(s in paste0("somite", c("I","II","III"))){
    df <- data.frame(var1=factor(c(rep("8",5), rep("18",4), rep("21",3), rep("25",2), "27"), levels = c("8", "18", "21", "25", "27")),
                     var2=factor(c("18","21","25","27","35", "21","25","27","35", "25","27","35", "27","35", "35"), levels = c("18", "21", "25", "27", "35")),
                     fc=as.numeric(somite[[s]][,2:16]),
                     fdr=somite[[s]]$FDR)
    df <- df[order(abs(df$fc)),]
    plots[[i]] <- ggplot(df, aes(var1, var2, fill = fc)) + geom_tile(color = ifelse(abs(df$fc) > log2(1.5), "black", "white"), size = ifelse(abs(df$fc) > log2(1.5), 0.5, 0)) + 
      scale_fill_gradient2(low = "steelblue", high = "indianred3", mid = "white", midpoint = 0, limit = rng, name="log2 FC") + 
      theme_minimal() + coord_fixed() + xlab(s) + ylab("") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text = element_text(face=2, size=10), axis.title.x = element_text(face=2, size=12))
    i <- i+1
  }
  p <- ggarrange(plotlist = plots, ncol = 1, nrow = 4, common.legend = TRUE, legend = "none")
  return(p)
  }
}


#### Somite trios
## Differential expression results
# MA plot
plotMAtriosAll <- function(contrast=3){
  tmp <- DEtriosAll[[as.numeric(contrast)]]
  colnames(tmp)[c(2,3,6)] <- c("log2FoldChange","baseMean","padj")
  tmp$baseMean <- 2^tmp$baseMean-1
  p <- ggmaplot(tmp, FDR = 0.05, fc=1.5, genenames = tmp$genes, size=1, legend = "top", top=30, 
                label.rectangle = TRUE, select.top.method = "padj", font.label = c("bold", 9), 
                main=names(DEtriosAll)[as.numeric(contrast)], ggtheme = theme_minimal() )
  return(p)
}

plotMAtriosStage <- function(contrast=1){
  tmp <- DEtriosStage[[as.numeric(contrast)]]
  colnames(tmp)[c(9,5,8)] <- c("log2FoldChange","baseMean","padj")
  tmp$baseMean <- 2^tmp$baseMean-1
  p <- ggmaplot(tmp, FDR = 0.05, fc=1.5, genenames = tmp$genes, size=1, legend = "top", top=30, 
                label.rectangle = TRUE, select.top.method = "padj", font.label = c("bold", 9), 
                main=names(DEtriosStage)[as.numeric(contrast)], ggtheme = theme_minimal() )
  return(p)
}

# DE results table
printDEtableTriosAll <- function(contrast=3){
  table <- DEtriosAll[[as.numeric(contrast)]][,-c(3:5)]
  table <- table[table$FDR < 0.05 & abs(table$logFC) > log2(1.5),]
  table$logFC <- round(table$logFC, 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  return(table)
}

printDEtableTriosStage <- function(contrast=1){
  table <- DEtriosStage[[as.numeric(contrast)]][,-c(5:7)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  table <- table[table$stageSpecific==1,]
  table <- table[,-c((ncol(table)-1):ncol(table))]
  table[,2] <- round(table[,2], 2)
  table[,3] <- round(table[,3], 2)
  table[,4] <- round(table[,4], 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  return(table)
}

# plot selected gene from DE table
boxplotExprTriosAll <- function(group_by=4, colour_by=3, contrast=3, selected=NULL){
  table <- DEtriosAll[[as.numeric(contrast)]][,-c(3:5)]
  table <- table[table$FDR < 0.05 & abs(table$logFC) > log2(1.5),]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 3){
      p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
    }
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 2){
      p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
    }
  }
  return(p)
}

boxplotExprTriosStage <- function(group_by=4, colour_by=3, contrast=1, selected=NULL){
  table <- DEtriosStage[[as.numeric(contrast)]][,-c(5:7)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  table <- table[table$stageSpecific==1,]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 3){
      p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
    }
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 2){
      p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
    }
  }
  return(p)
}

# GO enrichment results
printGOtableTrios <- function(level=1){
  table <- GOresultsTrios[GOresultsTrios$Fisher.classic.adj < 0.1,-6]
  colnames(table)[6] <- "FDR"
  table$FDR <- format(table$FDR, digits=4, width=4)
  row.names(table) <- table$GO.ID
  table$GO.ID <- NULL
  return(table)
}

# genes associated with selected GO term
getGOgenesTrios <- function(level=1, selected=NULL){
  table <- GOresultsTrios[GOresultsTrios$Fisher.classic.adj < 0.1,-6]
  colnames(table)[6] <- "FDR"
  table$FDR <- format(table$FDR, digits=4, width=4)
  row.names(table) <- table$GO.ID
  table$GO.ID <- NULL
  
  go <- row.names(table)[selected]
  table <- GOmapping[GOmapping$go_id==go,1]
  table <- as.character(normCounts[table,1])
  return(table)
}

boxplotExprTriosGO <- function(group_by=4, colour_by=3, gene=NULL){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 3){
      p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
    }
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
    if(group_by == 2){
      p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
    }
  }
  return(p)
}

#### Across development
## Differential expression results

# DE results table
printDEtableStageAll <- function(){
  table <- DEstageAll[,c(1,17,21,20)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  table$logCPM <- round(table$logCPM, 2)
  table$logFC.max <- round(table$logFC.max, 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  return(table)
}

printDEtableStageSomite <- function(contrast=1){
  table <- DEstageSomite[[as.numeric(contrast)]][,c(1,17,22,20,23)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  table <- table[table$somiteSpecific==1,]
  table$somiteSpecific <- NULL
  table[,2] <- round(table[,2], 2)
  table[,3] <- round(table[,3], 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  return(table)
}

# plot selected gene from DE table
boxplotExprStageAll <- function(group_by=4, colour_by=2, selected=NULL){
  table <- DEstageAll[,c(1,17,21,20)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
  }
  if(group_by == 2){
    fit <- lowess(d$cpm~d$stage)
    fit <- data.frame(x=fit$x, y=fit$y)
    p <- p + geom_line(data=fit, aes(x,y), alpha=0.5, lwd=1)
  }
  if(group_by == 4){
    fit <- lowess(d$cpm~d$somiteNumber)
    fit <- data.frame(x=fit$x, y=fit$y)
    p <- p + geom_line(data=fit, aes(x,y), alpha=0.5, lwd=1)
  }
  if(group_by == 2 & colour_by == 3){
    p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
  }
  if(group_by == 3 & colour_by == 2){
    p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
  }
  return(p)
}

boxplotExprStageSomite <- function(group_by=3, colour_by=2, contrast=1, selected=NULL){
  table <- DEstageSomite[[as.numeric(contrast)]][,c(1,17,22,20,23)]
  table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
  table <- table[table$somiteSpecific==1,]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "somiteNumber"))
  if(colour_by == 2){
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=col.stage) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
  }else{
    p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_manual(values=alpha(rep("orchid4",3),c(0.75, 0.5, 0.25))) + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + th
  }
  if(group_by == 2){
    fit <- lowess(d$cpm~d$stage)
    fit <- data.frame(x=fit$x, y=fit$y)
    p <- p + geom_line(data=fit, aes(x,y), alpha=0.5, lwd=1)
  }
  if(group_by == 4){
    fit <- lowess(d$cpm~d$somiteNumber)
    fit <- data.frame(x=fit$x, y=fit$y)
    p <- p + geom_line(data=fit, aes(x,y), alpha=0.5, lwd=1)
  }
  if(group_by == 2 & colour_by == 3){
    p <- p + geom_vline(xintercept = seq(1.5,5.5,1), lty=3, col="grey")
  }
  if(group_by == 3 & colour_by == 2){
    p <- p + geom_vline(xintercept = c(1.5,2.5), lty=3, col="grey")
  }
  return(p)
}













# GO enrichment results
printGOtableStage <- function(){
  table <- GOresultsStage
  table <- table[table$FDR < 0.1,-5]
  table$FDR <- format(table$FDR, digits=4, width=4)
  table$Term <- paste0(table$Term, " [", table$Ont, "]")
  table$Ont <- NULL
  return(table)
}

# genes associated with selected GO term
getGOgenesStage <- function(selected=NULL){
  table <- GOresultsStage
  table <- table[table$FDR < 0.1,-5]
  
  go <- row.names(table)[selected]
  table <- GOmapping[GOmapping$go_id==go,1]
  table <- as.character(allDEgenesStage[allDEgenesStage$id %in% table,1])
  return(table)
}

boxplotExprStageGO <- function(group_by=2, colour_by=3, gene=NULL){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=log10(cpm+1))) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ggtitle(gene) + theme(legend.position="none", axis.title = element_text(size = 15), axis.text = element_text(size=12))
  return(p)
}

## Hox correlated genes
printHoxCorrGenes <- function(){
  table <- hox[,c(8,7,3,5)]
  table$rho <- round(table$rho,2)
  table$FDR <- format(table$FDR, digits=4, width=4)
  return(table)
}

plotHoxCorrGenes <- function(idx=1){
  table <- hox[,c(8,7,3,5)]
  g1 <- row.names(normCounts[normCounts$gene == table[idx,1],])
  g2 <- row.names(normCounts[normCounts$gene == table[idx,2],])
  df <- data.frame(gene1 = as.numeric(log10(normCounts[g1,-1]+1)), gene2 = as.numeric(log10(normCounts[g2,-1]+1)), stage = meta$stage)
  ggplot(df, aes(gene1, gene2, colour=stage)) + geom_point(size=3) + scale_color_manual(values=levels(meta$col)) + xlab(normCounts[g1,1]) + ylab(normCounts[g2,1]) + theme(axis.title = element_text(size=15))
}

printGOtableHox <- function(){
  table <- GOresultsHox
  table <- table[table$FDR < 0.1,-5]
  table$FDR <- format(table$FDR, digits=4, width=4)
  table$Term <- paste0(table$Term, " [", table$Ont, "]")
  table$Ont <- NULL
  return(table)
}

getGOgenesHox <- function(selected=NULL){
  table <- GOresultsHox
  table <- table[table$FDR < 0.1,-5]
  
  go <- row.names(table)[selected]
  table <- GOmapping[GOmapping$go_id==go,1]
  table <- as.character(allDEgenesHox[allDEgenesHox$id %in% table,1])
  return(table)
}

boxplotExprHoxGO <- function(gene=NULL){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,2], y=log10(cpm+1))) + geom_boxplot(aes(fill=d[,3])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[2]) + ggtitle(gene) + theme(legend.position="none", axis.title = element_text(size = 15), axis.text = element_text(size=12))
  return(p)
}

##################################################
### prepare environment
# library(biomaRt)
# dir <- "/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/"
# ## metadata
# meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), header = TRUE, stringsAsFactors = FALSE)
# meta <- meta[meta$use==1,]
# meta$somite <- factor(meta$somite, c("SIII","SII","SI"))
# meta$stage <- as.factor(meta$stage)
# meta$date <- as.factor(meta$date)
# meta$group <- as.factor(meta$group)
# meta$somiteNumber <- factor(meta$somiteNumber, levels=c(6:8,16:21,23:27,33:35))
# 
# ## expression estimates
# normCounts <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_batchCorrected_14PCs.tsv"))
# identical(meta$sample, colnames(normCounts)[-1])
# 
# ## gene - GO mappings
# ensembl <- useMart(host='apr2019.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="mmusculus_gene_ensembl") # v96
# GOmapping <- getBM(attributes=c('ensembl_gene_id', 'go_id'), filters = 'ensembl_gene_id', values = row.names(normCounts), mart = ensembl)
# GOmapping <- GOmapping[GOmapping$go_id != "",]
# 
# ## somite trios
# DEtriosAll <- list()
# for(contrast in paste0("somite", c("IvsII", "IIvsIII", "IvsIII"))){
#   DEtriosAll[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/03.3_DEresults_somiteTrios_", contrast, "_pca.tsv"), stringsAsFactors = FALSE)
# }
# # only keep in per-stage those that are not in the average test
# tmp <- c(row.names(DEtriosAll[["somiteIvsII"]][DEtriosAll[["somiteIvsII"]]$FDR<0.05 & abs(DEtriosAll[["somiteIvsII"]]$logFC) > log2(1.5),]),
#          row.names(DEtriosAll[["somiteIIvsIII"]][DEtriosAll[["somiteIIvsIII"]]$FDR<0.05 & abs(DEtriosAll[["somiteIIvsIII"]]$logFC) > log2(1.5),]),
#          row.names(DEtriosAll[["somiteIvsIII"]][DEtriosAll[["somiteIvsIII"]]$FDR<0.05 & abs(DEtriosAll[["somiteIvsIII"]]$logFC) > log2(1.5),]) )
# tmp <- unique(tmp)
# DEtriosStage <- list()
# for(contrast in paste0("stage", c(8,18,21,25,27,35))){
#   DEtriosStage[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/03.3_DEresults_somiteTrios_", contrast, "_pca.tsv"), stringsAsFactors = FALSE)
#   DEtriosStage[[contrast]]$logFC.max <- sapply(1:nrow(DEtriosStage[[contrast]]), function(x) DEtriosStage[[contrast]][x,which.max(abs(DEtriosStage[[contrast]][x,2:4]))+1])
#   DEtriosStage[[contrast]]$stageSpecific <- ifelse(DEtriosStage[[contrast]]$FDR < 0.05 & abs(DEtriosStage[[contrast]]$logFC.max) > log2(1.5) & !(row.names(DEtriosStage[[contrast]]) %in% tmp), 1, 0)
# }
# GOresultsTrios <- read.table(paste0(dir, "RNA-seq/results/03.5_DEgenes_somiteTrios_all_GOresults.tsv"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
# 
# ## stages
# DEstageAll <- read.table(paste0(dir, "RNA-seq/results/03.3_DEresults_stage_all_pca.tsv"), stringsAsFactors = FALSE)
# DEstageAll$logFC.max <- sapply(1:nrow(DEstageAll), function(x) DEstageAll[x,which.max(abs(DEstageAll[x,2:16]))+1])
# # keep as somite-specific only those not in average test
# tmp <- row.names(DEstageAll[DEstageAll$FDR<0.05 & abs(DEstageAll$logFC.max) > log2(1.5),])
# DEstageSomite <- list()
# for(contrast in paste0("somite", c("I", "II", "III"))){
#   DEstageSomite[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/03.3_DEresults_stage_",contrast,"_pca.tsv"), stringsAsFactors = FALSE)
#   DEstageSomite[[contrast]]$logFC.max <- sapply(1:nrow(DEstageSomite[[contrast]]), function(x) DEstageSomite[[contrast]][x,which.max(abs(DEstageSomite[[contrast]][x,2:16]))+1])
#   DEstageSomite[[contrast]]$somiteSpecific <- ifelse(DEstageSomite[[contrast]]$FDR<0.05 & abs(DEstageSomite[[contrast]]$logFC.max) > log2(1.5) & !(row.names(DEstageSomite[[contrast]]) %in% tmp), 1, 0)
# }
# GOresultsStage <- read.table(paste0(dir, "RNA-seq/results/03.5_DEgenes_stages_all_GOresults.tsv"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
# 
# # hox <- read.table("OneDrive/OneDrive - Cancer Research UK, Cambridge Institute/SOMITES/RNAseq/HOX/HoxCorrelatedGenes.tab", stringsAsFactors = FALSE)
# # colnames(hox)[7] <- "geneName2"
# # hox$geneName1 <- as.character(normCounts[hox$gene1,1])
# # allDEgenesHox <- read.table("OneDrive/OneDrive - Cancer Research UK, Cambridge Institute/SOMITES/RNAseq/HOX/corrGenesGOtesting.tab")
# # colnames(allDEgenesHox) <- c("name", "id")
# #
# # GOresultsHox <- read.table("OneDrive/OneDrive - Cancer Research UK, Cambridge Institute/SOMITES/RNAseq/HOX/GOenrichment_HoxCorrelatedGenes.tab", sep="\t")
# 
# ## save
# save(meta, normCounts, GOmapping, DEtriosAll, DEtriosStage, GOresultsTrios, DEstageSomite, DEstageAll, GOresultsStage, file=paste0(dir, "shinyApps/RNA-seq/data/data.RData"))







