# helper.R
load("data/data.RData")

retrieveGeneExpr <- function(normCounts, meta, columns=c("stage", "somite", "date"), gene){
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
boxplotExpr <- function(group_by=2, colour_by=3, gene="Hoxa1"){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + theme( axis.title = element_text(size = 15), axis.text = element_text(size=12))
  return(p)
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
  table <- table[table$FDR < 0.05,]
  table$logFC <- round(table$logFC, 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  if(sum(duplicated(table$genes))>0){
    table$genes <- as.character(table$genes)
    idx <- which(duplicated(table$genes))
    g <- table[idx,]$genes
    idx <- which(table$genes == g)
    count <- 0
    for(dup in g){
      count <- count+1;
      j <- 1
      for(i in idx[(2*count-1):(2*count)]){
        table[i,]$genes <- paste(dup, j, sep=".")
        j <- j+1
      }
    }
  }
  row.names(table) <- table$genes
  table$genes <- NULL
  return(table)
}

printDEtableTriosStage <- function(contrast=1){
  table <- DEtriosStage[[as.numeric(contrast)]][,-c(5:7,9)]
  table <- table[table$FDR < 0.05,]
  table <- table[table$stageSpecific==1,]
  table$stageSpecific <- NULL
  table[,2] <- round(table[,2], 2)
  table[,3] <- round(table[,3], 2)
  table[,4] <- round(table[,4], 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  if(sum(duplicated(table$genes))>0){
    table$genes <- as.character(table$genes)
    idx <- which(duplicated(table$genes))
    g <- table[idx,]$genes
    idx <- which(table$genes == g)
    count <- 0
    for(dup in g){
      count <- count+1;
      j <- 1
      for(i in idx[(2*count-1):(2*count)]){
        table[i,]$genes <- paste(dup, j, sep=".")
        j <- j+1
      }
    }
  }
  row.names(table) <- table$genes
  table$genes <- NULL
  return(table)
}

# plot selected gene from DE table
boxplotExprTriosAll <- function(group_by=2, colour_by=3, contrast=3, selected=NULL){
  table <- DEtriosAll[[as.numeric(contrast)]][,-c(4:5)]
  table <- table[table$FDR < 0.05,]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + theme(axis.title = element_text(size = 15), axis.text = element_text(size=12))
  return(p)
}

boxplotExprTriosStage <- function(group_by=2, colour_by=3, contrast=1, selected=NULL){
  table <- DEtriosStage[[as.numeric(contrast)]][,-c(5:7)]
  table <- table[table$FDR < 0.1,]
  table <- table[table$stageSpecific==1,]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)])
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

boxplotExprTriosGO <- function(group_by=2, colour_by=3, gene=NULL){
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + theme(axis.title = element_text(size = 15), axis.text = element_text(size=12))
  return(p)
}

#### Across development
## Differential expression results

# DE results table
printDEtableStageAll <- function(){
  table <- DEstageAll[,c(1,17,21,20)]
  table <- table[table$FDR < 0.05,]
  table <- table[abs(table$logFC.max) > log2(1.5),]
  table$logCPM <- round(table$logCPM, 2)
  table$logFC.max <- round(table$logFC.max, 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  if(sum(duplicated(table$genes))>0){
    table$genes <- as.character(table$genes)
    idx <- which(duplicated(table$genes))
    g <- table[idx,]$genes
    idx <- which(table$genes %in% g)
    count <- 0
    for(dup in g){
      count <- count+1;
      j <- 1
      for(i in idx[(2*count-1):(2*count)]){
        table[i,]$genes <- paste(dup, j, sep=".")
        j <- j+1
      }
    }
  }
  row.names(table) <- table$genes
  table$genes <- NULL
  return(table)
}

printDEtableStageSomite <- function(contrast=1){
  table <- DEstageSomite[[as.numeric(contrast)]][,c(1,17,22,20,23)]
  table <- table[table$FDR < 0.05,]
  table <- table[table$somiteSpecific==1,]
  table$somiteSpecific <- NULL
  table[,2] <- round(table[,2], 2)
  table[,3] <- round(table[,3], 2)
  table$FDR <- format(table$FDR, digits=4, width=5)
  if(sum(duplicated(table$genes))>0){
    table$genes <- as.character(table$genes)
    idx <- which(duplicated(table$genes))
    g <- table[idx,]$genes
    idx <- which(table$genes == g)
    count <- 0
    for(dup in g){
      count <- count+1;
      j <- 1
      for(i in idx[(2*count-1):(2*count)]){
        table[i,]$genes <- paste(dup, j, sep=".")
        j <- j+1
      }
    }
  }
  row.names(table) <- table$genes
  table$genes <- NULL
  return(table)
}

# printDEtableStagePairwise <- function(contrast=1){
#   col <- as.numeric(contrast)+1
#   table <- DEstage[,c(1,17,col,20)]
#   table <- table[table$FDR < 0.1,]
#   table$logCPM <- round(table$logCPM, 2)
#   table[,3] <- round(table[,3], 2)
#   table$FDR <- format(table$FDR, digits=4, width=5)
#   table <- table[order(abs(table[,3]), decreasing = TRUE),]
#   if(sum(duplicated(table$genes))>0){
#     table$genes <- as.character(table$genes)
#     idx <- which(duplicated(table$genes))
#     g <- table[idx,]$genes
#     idx <- which(table$genes == g)
#     count <- 0
#     for(dup in g){
#       count <- count+1;
#       j <- 1
#       for(i in idx[(2*count-1):(2*count)]){
#         table[i,]$genes <- paste(dup, j, sep=".")
#         j <- j+1
#       }
#     }
#   }
#   row.names(table) <- table$genes
#   table$genes <- NULL
#   return(table)
# }

# plot selected gene from DE table
boxplotExprStageAll <- function(group_by=2, colour_by=2, selected=NULL){
  table <- DEstageAll[,c(1,17,21,20)]
  table <- table[table$FDR < 0.05,]
  table <- table[abs(table$logFC.max) > log2(1.5),]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  fit <- lowess(d$cpm~d$stage)
  fit <- data.frame(x=fit$x, y=fit$y)
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)]) + theme(axis.title = element_text(size = 15), axis.text = element_text(size=12))
  if(group_by==2){ p <- p + geom_line(data=fit, aes(x,y), col="mediumorchid4", lwd=1) }
  return(p)
}

boxplotExprStagePairwise <- function(group_by=3, colour_by=2, contrast=1, selected=NULL){
  table <- DEstageSomite[[as.numeric(contrast)]][,c(1,17,22,20,23)]
  table <- table[table$FDR < 0.05,]
  table <- table[table$somiteSpecific==1,]
  gene <- as.character(table[selected,1])
  
  d <- retrieveGeneExpr(normCounts, meta, gene=gene, columns=c("stage", "somite", "date"))
  fit <- lowess(d$cpm~d$stage)
  fit <- data.frame(x=fit$x, y=fit$y)
  p <- ggplot(d, aes(x=d[,as.numeric(group_by)], y=cpm)) + geom_boxplot(aes(fill=d[,as.numeric(colour_by)])) + scale_fill_brewer(palette = "Purples") + xlab(colnames(d)[as.numeric(group_by)]) + ylab("log2 CPM") + ggtitle(gene) + labs(fill=colnames(d)[as.numeric(colour_by)])
  if(group_by==2){ p <- p + geom_line(data=fit, aes(x,y), col="mediumorchid4", lwd=1) }
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
# dir <- "/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/"
# ## metadata
# meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), header = TRUE, stringsAsFactors = FALSE)
# meta <- meta[meta$use==1,]
# meta$somite <- as.factor(meta$somite)
# meta$stage <- as.factor(meta$stage)
# meta$date <- as.factor(meta$date)
# meta$group <- as.factor(meta$group)
# 
# ## expression estimates
# normCounts <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"))
# ## remove technical variation
# pcs <- read.table(paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"))
# tmp <- removeBatchEffect(normCounts[,-1], design = model.matrix(~0+group, meta), covariates = pcs[,1:14])
# normCounts <- cbind(normCounts[match(row.names(tmp), row.names(normCounts)),1], as.data.frame(tmp))
# colnames(normCounts)[1] <- "gene"
# rm(pcs); rm(tmp)
# 
# identical(meta$sample, colnames(normCounts)[-1])
# 
# ## gene - GO mappings
# ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# GOmapping <- getBM(attributes=c('ensembl_gene_id', 'go_id'), filters = 'ensembl_gene_id', values = row.names(normCounts), mart = ensembl)
# GOmapping <- GOmapping[GOmapping$go_id != "",]
# 
# ## somite trios
# DEtriosAll <- list()
# for(contrast in paste0("somite", c("IvsII", "IIvsIII", "IvsIII"))){
#   DEtriosAll[[contrast]] <- read.table(paste0(dir, "RNA-seq/results/03.3_DEresults_somiteTrios_", contrast, "_pca.tsv"), stringsAsFactors = FALSE)
# }
# # only keep in per-stage those that are not in the average test
# tmp <- c(row.names(DEtriosAll[["somiteIvsII"]][DEtriosAll[["somiteIvsII"]]$FDR<0.05 & abs(DEtriosAll[["somiteIvsII"]]$logFC) > log2(1.5),]), row.names(DEtriosAll[["somiteIIvsIII"]][DEtriosAll[["somiteIIvsIII"]]$FDR<0.05 & abs(DEtriosAll[["somiteIIvsIII"]]$logFC) > log2(1.5),]),
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
# keep as somite-specific only those not in average test
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






