---
title: "Removal of technical effects from RNA-seq data"
date: "20 May 2019"
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: true
    toc_float: 
      collapsed: false
---



We have processed and QCed RNA-seq data from mouse somites. All but one sample are of good enough quality for downstream analyses. We ignore the `PSM` samples, since these were to assess any problems with particular batches of the data, but we don't have any serious failures. From now on, we focus on 76 good-quality somite samples.


```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$use==1,]

data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"))
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"))

y <- readRDS(paste0(dir, "RNA-seq/results/01_edgeRobject.Rds"))
```

We observed strong technical effects relating to the date of sample collection. Since the experimental design is confounded with collection date, we cannot simply regress it out on the differential expression analysis. Thus, we need to somehow identify and remove these technical variation.

### SVA

One approach is to use surrogate variable analysis. We can run `svaseq` on the normalised counts, providing no *known* batch effects. The idea is that the surrogate variables will capture the technical batch from collection date and we can use these SVs to regress it out.

The algorithm identifies 11 surrogate variables. The first one correlates pretty well with the date of sample collection.


```r
## design matrix
meta$group <- as.factor(meta$group)
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage", levels(meta$group))

## SVA
mod0 = model.matrix(~1, data=meta)
svobj = svaseq(cpm(y), design, mod0)
```

```
## Number of significant surrogate variables is:  11 
## Iteration (out of 5 ):1  2  3  4  5
```

```r
tmp <- svobj$sv
tmp <- cbind(tmp, meta)
tmp$date <- factor(as.character(tmp$date), levels = c("06.04.18", "05.04.18", "12.04.18", "20.05.18", "26.04.18", "24.05.18", "01.06.18", "11.06.18", "27.04.18", "14.06.18"))
tmp$stage <- as.factor(tmp$stage)
tmp$somite <- as.factor(tmp$somite)

par(mfrow=c(1,3))
for(i in 1:ncol(svobj$sv)){
  plot(as.numeric(tmp$date), tmp[,i], pch=16, col=tmp$date, axes=FALSE, xlab="", ylab=paste0("SV",i), cex=2); box(bty="l"); axis(1, at=1:10, labels = levels(tmp$date), las=2)
  plot(as.numeric(tmp$stage), tmp[,i], pch=16, col=tmp$stage, axes=FALSE, xlab="", ylab=paste0("SV",i), cex=2); box(bty="l"); axis(1, at=1:6, labels = levels(tmp$stage), las=2)
  plot(as.numeric(tmp$somite), tmp[,i], pch=16, col=tmp$somite, axes=FALSE, xlab="", ylab=paste0("SV",i), cex=2); box(bty="l"); axis(1, at=1:3, labels = levels(tmp$somite), las=2)
}
```

![](02_removeTechnicalVar_files/figure-html/sva-1.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-2.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-3.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-4.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-5.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-6.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-7.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-8.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-9.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-10.png)<!-- -->![](02_removeTechnicalVar_files/figure-html/sva-11.png)<!-- -->

With so many SVs it is a bit difficult to know whether some are capturing biological variation not explicit in the model, which could be unfortunate.

To investigate this, we look at whether any genes are strongly correlated with the SVs. Ideally, most genes won't be.


```r
gene.lists <- list(); j=1
par(mfrow=c(3,4))
for(i in 1:ncol(svobj$sv)){
  cors <- sapply(1:nrow(dataNorm), function(x) cor(svobj$sv[,i], as.numeric(dataNorm[x,-1])))
  g <- as.character(y$genes[which(abs(cors)>0.75),])
  plot(density(cors), main=paste0("SV",i))
  mtext(side=3, line=-1, text = paste(length(g), "genes with abs(r)>0.75"), cex=0.75, adj=0)
  if(length(g)>40){
    gene.lists[[j]] <- g
    j <- j+1
  }
}
```

![](02_removeTechnicalVar_files/figure-html/corrGenesSVs-1.png)<!-- -->

The first 3 SVs are correlated to some genes. We check whether these genes are involved in articular functions, or are they just random sets. We run a gene ontology enrichment on the three gene lists.


```r
universe <- row.names(dataNorm)

go <- list()
go.test <- list()
go.res <- list()
for(i in 1:length(gene.lists)){
  sel <- row.names(dataNorm[dataNorm$gene %in% gene.lists[[i]],])
  all <- as.factor(as.numeric(universe %in% sel))
  names(all) <- universe
  go[[i]] <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
  go.test[[i]] <- runTest(go[[i]], algorithm = "classic", statistic = "Fisher" )
  go.res[[i]] <- GenTable(go[[i]], Fisher.classic = go.test[[i]], topNodes = length(score(go.test[[i]])))
  go.res[[i]]$Fisher.classic.adj <- p.adjust(go.res[[i]]$Fisher.classic, "fdr")
}

# go.res[[1]] # adhesion; vasculature; angiogenesis; locomotion; signaling
# go.res[[2]] # translation; metabolic processes; 
# go.res[[3]] # gene expression; (translation, metabolism)
```

- The first SV, correlated with 116 genes that are enriched for terms related to vasculature and angiogenesis. This suggests that we might be picking up the amount of *contamination* from blood and endothelial cells. If this is the case, it is ok to regress this out.
- The second SV, correlated with 666 genes that are enriched for terms related to translation, transcription and macromolecule metabolic processes.
- The third SV, correlated with 47 genes that are not significantly enriched for any terms.

Thus, it seems that the SVs are not capturing any obvious biological processes that we are interested in.

We can check if the batch effect has been successfully removed by regressing out the SVs and redoing the PCA on the top 1000 most variable genes.

Samples cluster very clearly by stage, and now we can also see grouping by somite. The grouping by date is no longer evident.


```r
dataNorm.sva <- removeBatchEffect(dataNorm[,-1], design = model.matrix(~0+group, meta), covariates = svobj$sv)

vars <- rowVars(dataNorm.sva)
names(vars) <- row.names(dataNorm.sva)
pca <- prcomp(t(dataNorm.sva[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$stage))) + labs(col="stage") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$somite))) + labs(col="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[3]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$date))) + labs(col="batch") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
multiplot(plotlist = plots, cols=3)
```

![](02_removeTechnicalVar_files/figure-html/plotPCA-1.png)<!-- -->


### PCA

An alternative approach is to use PCA on the residuals after fitting the model including our biological effects and remove the amount of variation that is above that expected by chance. This will only preserve what is explicitly modelled in our design, that is the *stage* and *somite* of each sample; any other effects, technical *and* biological, will be removed.

We fit a model including the design considering stage and somite, and apply PCA on the residuals of the fit. We then use the `parallelPCA` function from `scran` to estimate the number of PCs to keep. This test is based on permuting the data to establish how much variation is explained under a random model; PCs that capture more variation than the random model are kept.


```r
## insead of using SVA, use PCA to identify the variation on the residuals after fitting the design of interest
fit <- lmFit(dataNorm[,-1], design)
res <- residuals(fit, dataNorm[,-1])

# pca
pcs <- prcomp(t(res))
```

For this dataset we retain the first 14 PCs.

Again, we can regress out these PCs and rerun the PCA to check the batch effect. The separation by stage is striking, and within each stage there is pretty clear separation by somite also. Clustering by batch is no longer observed.


```r
# plugging 'dataNorm' into scran::parallelPCA() returns 14
dataNorm.pca <- removeBatchEffect(dataNorm[,-1], design = model.matrix(~0+group, meta), covariates = pcs$x[,1:14])

vars <- rowVars(dataNorm.pca)
names(vars) <- row.names(dataNorm.pca)
pca <- prcomp(t(dataNorm.pca[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$stage))) + labs(col="stage") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$somite))) + labs(col="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[3]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$date))) + labs(col="batch") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
multiplot(plotlist = plots, cols=3)
```

![](02_removeTechnicalVar_files/figure-html/plotPCA2-1.png)<!-- -->

So this approach cleans out the data nicely and allows us already to reconstruct a trajectory across the AP axis, from the earliest to the latest somites:


```r
ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.numeric(meta[meta$use==1,]$somiteNumber))) + labs(col="somite number") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained")) + scale_color_gradientn(colours=brewer.pal(n=9, "Blues")[-c(1:2)])
```

![](02_removeTechnicalVar_files/figure-html/somiteNumber-1.png)<!-- -->

We can now use either the SVs or PC from the residuals as covariates to regress out when performing differential expression, to ensure that the technical effects are properly controlled for.


```r
write.table(svobj$sv, paste0(dir, "RNA-seq/results/02_svs.tab"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(pcs$x, paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
```



```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS: /datastore/sw/gentoo2/usr/lib64/blas/reference/libblas.so.3.7.0
## LAPACK: /datastore/sw/gentoo2/usr/lib64/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] org.Mm.eg.db_3.7.0   topGO_2.34.0         SparseM_1.77        
##  [4] GO.db_3.7.0          AnnotationDbi_1.44.0 IRanges_2.16.0      
##  [7] S4Vectors_0.20.1     Biobase_2.42.0       graph_1.60.0        
## [10] BiocGenerics_0.28.0  sva_3.30.1           BiocParallel_1.16.6 
## [13] genefilter_1.64.0    mgcv_1.8-24          nlme_3.1-137        
## [16] RColorBrewer_1.1-2   ggplot2_3.1.1        edgeR_3.24.3        
## [19] limma_3.38.3        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1         locfit_1.5-9.1     lattice_0.20-35   
##  [4] assertthat_0.2.1   digest_0.6.18      R6_2.4.0          
##  [7] plyr_1.8.4         RSQLite_2.1.1      evaluate_0.13     
## [10] pillar_1.3.1       rlang_0.3.4        lazyeval_0.2.2    
## [13] annotate_1.60.1    blob_1.1.1         Matrix_1.2-14     
## [16] rmarkdown_1.12     labeling_0.3       splines_3.5.1     
## [19] stringr_1.4.0      RCurl_1.95-4.12    bit_1.1-14        
## [22] munsell_0.5.0      compiler_3.5.1     xfun_0.6          
## [25] pkgconfig_2.0.2    htmltools_0.3.6    tidyselect_0.2.5  
## [28] tibble_2.1.1       matrixStats_0.54.0 XML_3.98-1.15     
## [31] crayon_1.3.4       dplyr_0.8.0.1      withr_2.1.2       
## [34] bitops_1.0-6       xtable_1.8-3       gtable_0.3.0      
## [37] DBI_1.0.0          magrittr_1.5       scales_1.0.0      
## [40] stringi_1.4.3      tools_3.5.1        bit64_0.9-7       
## [43] glue_1.3.1         purrr_0.3.2        survival_2.42-6   
## [46] yaml_2.2.0         colorspace_1.4-1   memoise_1.1.0     
## [49] knitr_1.22
```

