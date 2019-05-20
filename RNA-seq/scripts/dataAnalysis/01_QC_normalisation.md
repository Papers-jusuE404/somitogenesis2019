---
title: "RNA-seq data QC and normalisation"
date: "17 May 2019"
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



We have sequenced 80 RNA-seq libraries from individual mouse somites, generally matched to the ATAC-seq samples; in some cases, library construction was only successful for one technique giving rise to some singletons, but mostly, samples contain both data types. Additionally, we included 12 samples of *adjacent* tissue as controls.

Samples represent the three most recently segmented somites per embryo, across six different stages; stages are named based on the total number of somites in the embryo. We have 3 to 6 biological replicates per stage.


```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
table(meta$somite, meta$stage)[-1,]
```

```
##       
##        8 18 21 25 27 35
##   SI   4  5  4  4  4  6
##   SII  4  5  3  3  4  5
##   SIII 4  5  4  4  3  6
```

We have mapped the data to the mouse reference genome (`mm10`) and quantified expression by counting the number of aligned fragments per gene (using Ensembl annotation, version 96). 

### Quality control

Samples were sequenced at a median depth of nearly 17.5 million fragments; 75% of the data have library sizes between 14.9 and 23.2 million fragments.


```r
## data
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), stringsAsFactors = FALSE)

## match column order in count matrix and metadata annotation
meta <- meta[match(colnames(data)[-1], meta$sample),]
stopifnot(identical(colnames(data)[-1], meta$sample))

## mapping statistics
mapping.stats <- read.table(paste0(dir, "RNA-seq/data/mappingStatistics.tsv"), stringsAsFactors = FALSE)

stopifnot(identical(colnames(data)[-1], row.names(mapping.stats)))
mapping.stats$uniquelyMapped <- colSums(data[,-1])
mapping.stats$libSize <- rowSums(mapping.stats)
# summary(mapping.stats$libSize)

ggplot(mapping.stats, aes(1, libSize/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("library size (million fragments)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_normalisation_files/figure-html/depth-1.png)<!-- -->

And around 12.8 million of these fragments are uniquely mapped and within exons (interquartile range = 10.5 - 16.4 million)


```r
# summary(mapping.stats$uniquelyMapped)
ggplot(mapping.stats, aes(1, uniquelyMapped/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("total fragments uniquely mapped to exons (millions)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_normalisation_files/figure-html/inExons-1.png)<!-- -->

One sample, e17_SII.2, has a very small library with only 8.843\times 10^{4} total fragments and will be removed from downstream analyses. The next smallest library has 5.3 million fragments which might or might not be enough for a representative transcriptome.

In general, most samples show good mapping statistics with low proportion of unmapped and multimapped fragments, and of fragments mapped to non-exonic regions. There are a few outliers that show larger numbers, but since these have large libraries (over 20 million), the proportion of uniquely mapped fragments is still large.


```r
plots <- list()
plots[[1]] <- ggplot(mapping.stats, aes(1, N_unmapped/libSize*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_unmapped/libSize*100, col=mapping.stats$libSize/1e6), width = 0.05) + ylab("% unmapped") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[2]] <- ggplot(mapping.stats, aes(1, N_multimapping/libSize*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_multimapping/libSize*100, col=mapping.stats$libSize/1e6), width = 0.05) + ylab("% multimapped") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[3]] <- ggplot(mapping.stats, aes(1, N_noFeature/libSize*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_noFeature/libSize*100, col=mapping.stats$libSize/1e6), width = 0.05) + ylab("% outside exons") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plots[[4]] <- ggplot(mapping.stats, aes(1, N_ambiguous/libSize*100)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,N_ambiguous/libSize*100, col=mapping.stats$libSize/1e6), width = 0.05) + ylab("% ambiguous") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

multiplot(plotlist = plots, cols=2)
```

![](01_QC_normalisation_files/figure-html/qc-1.png)<!-- -->

The two samples with the smallest libraries have good stats, with only ~7% of fragments unmapped, multimapped, outside exons and ambiguous. That leaves around 70% of fragments uniquely mapped and thus suggests the libraries are of good quality, but sequenced to insufficient depth.

Samples generally have around 22 thousand genes expressed, with a tight IQR between 21.1 and 23 thousand. The sample with the very small library is a clear outlier and needs to be removed. The second smallest also deviates a bit from the rest of the dataset, but has 18,565 detected genes, compared to 19,549 of the next sample.


```r
mapping.stats$nGenes <- apply(data[,-1], 2, function(x) sum(x>0))
# summary(mapping.stats$nGenes)

ggplot(mapping.stats, aes(1, nGenes/1e3)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,nGenes/1e3, col=mapping.stats$libSize/1e6), width = 0.02) + ylab("total genes detected (thousands)") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](01_QC_normalisation_files/figure-html/numGenes-1.png)<!-- -->

Given that the second smallest library shows a similar number of detected genes and good mapping statistics (reflecting good sample quality), it is likely that its data is useful, specially for genes expressed at moderate to high levels. Therefore, we will only remove the sample that almost failed sequencing and proceed with all others.


```r
bad.qual <- row.names(mapping.stats[which.min(mapping.stats$libSize),])
data <- data[,-which(colnames(data) == bad.qual)]
meta$QC <- ifelse(meta$sample == bad.qual, 0, 1)
```

### Normalisation

The first thing to account for is library size. Before estimating size factors, we filter very lowly expressed genes, since we don't have enough statistical power to compare their expression across conditions, and can make the size factor estimation unstable.

We retain around 40% of the annotated genes (46% of detected genes).


```r
## our factor of interest is the combination of stage and somite. 
# create a variable for it
meta$group <- paste(meta$stage, meta$somite, sep=".")

## create an edgeR object
y <- DGEList(counts=data[,-1], samples = meta[meta$QC==1,], genes = data[,1], group = meta[meta$QC==1,]$group)

## filter low abundance genes
means <- aveLogCPM(y)
keep <- filterByExpr(y)

plot(density(means[keep]), lwd=2, xlab="average log counts-per-million", main="", bty="l", ylim=c(0,1))
lines(density(means[-keep]), lty=2, lwd=2)
legend("topright", legend = c("kept", "filtered"), lty=c(1,2), lwd=2)
```

![](01_QC_normalisation_files/figure-html/filter-1.png)<!-- -->

```r
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   33675   21861
```

```r
y <- y[keep, , keep.lib.sizes=FALSE]
```

With this reduced dataset of 21861 genes we calculate size factors to normalise the data. Size factors normalise for both library size and composition biases.

The method successfully standardises the libraries. The sample with second smallest library size and slightly fewer detected genes is highlighted in red; the normalisation seems to successfully bring it to the overall dataset distribution, so it should be ok to keep.


```r
## normalisation
y <- calcNormFactors(y)

dataNorm <- cpm(y, log=TRUE, prior.count = 1)

col <- c(rep("black",18), "red", rep("black",69)) # highlight the sample with small library size
par(mfrow=c(1,2))
boxplot(log2(data[keep,-1]+1), main="RAW", ylab=expression('log'[2]*' counts + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
boxplot(dataNorm, main="NORMALISED", ylab=expression('log'[2]*' counts-per-million + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
```

![](01_QC_normalisation_files/figure-html/norm-1.png)<!-- -->

### Exploratory analysis

With the normalised data, we can explore the general features of the dataset. Firstly, we compute a PCA on the thousand most variable genes to check if we can detect our biological effect of interest and any potential confounders.

The first two PCs separate samples by their stage, suggesting there are fairly big differences in their transcriptomes. In contrast, there is no obvious grouping by somite, suggesting differences are much more subtle and perhaps captured in later PCs; but the PSM samples separate fairly well.

There is also evident grouping based on the date of collection. This is a technical batch effect that we need to account for.


```r
vars <- rowVars(dataNorm)
names(vars) <- row.names(dataNorm)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$stage))) + labs(col="stage") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$somite))) + labs(col="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
plots[[3]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$QC==1,]$date))) + labs(col="batch") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
multiplot(plotlist = plots, cols=3)
```

![](01_QC_normalisation_files/figure-html/pca-1.png)<!-- -->

Since the PSM samples are grouping better by their identity than their technical batch, it suggests that there were no big issues with any of the batches. That is, there is no anomalous behaviour of a group of samples based on their batch. Thus, we can remove the PSM samples and focus the analysis on the somite samples only.

### Somite dataset only

We remove all the PSM samples.


```r
psm <- meta[meta$somite=="PSM",]$sample
data <- data[,-which(colnames(data) %in% psm)]

meta$use <- ifelse(meta$sample %in% psm, 0, ifelse(meta$QC == 0, 0, 1))
```

We recalculate the genes to filter out, with similar numbers as before.


```r
y <- DGEList(counts=data[,-1], samples = meta[meta$use==1,], genes = data[,1], group = meta[meta$use==1,]$group)

## filter low abundance genes
means <- aveLogCPM(y)
keep <- filterByExpr(y)

# plot(density(means[keep]), lwd=2, xlab="average log counts-per-million", main="", bty="l", ylim=c(0,1))
# lines(density(means[-keep]), lty=2, lwd=2)
# legend("topright", legend = c("kept", "filtered"), lty=c(1,2), lwd=2)

summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   35474   20062
```

```r
y <- y[keep, , keep.lib.sizes=FALSE]
```
 
And the normalisation is again successful.
 

```r
## normalisation
y <- calcNormFactors(y)

dataNorm <- cpm(y, log=TRUE, prior.count = 1)

col <- c(rep("black",18), "red", rep("black",69)) # highlight the sample with small library size
par(mfrow=c(1,2))
boxplot(log2(data[keep,-1]+1), main="RAW", ylab=expression('log'[2]*' counts + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
boxplot(dataNorm, main="NORMALISED", ylab=expression('log'[2]*' counts-per-million + 1'), axes=FALSE, border=col); box(bty="l"); axis(2, las=2)
```

![](01_QC_normalisation_files/figure-html/norm2-1.png)<!-- -->

In a PCA we observe a similar separation by stage within the first two PCs, and clear grouping by collection batch.


```r
vars <- rowVars(dataNorm)
names(vars) <- row.names(dataNorm)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

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

![](01_QC_normalisation_files/figure-html/pca2-1.png)<!-- -->

The first PC is dominated by collection date, whereas stage gets pushed to PC2. Thus, we have strong batch effects from sample processing that need to be corrected.

However, the experimental design is confounded with batch. We have ten collection dates, and the six stages split into two blocks, each spanning five dates. Thus, blocking through the design matrix will not be possible.

For now, let's save the results.


```r
## metadata with the additional info of which samples are used in downstream analyses
write.table(meta, paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), quote = FALSE, sep="\t", row.names = FALSE)

## save the normalised counts, but first add gene names
dataNorm <- cbind(data[match(row.names(dataNorm), row.names(data)),1], as.data.frame(dataNorm))
colnames(dataNorm)[1] <- "gene"
write.table(dataNorm, paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"), quote = FALSE, sep="\t")

## also save edgeR object
saveRDS(y, paste0(dir, "RNA-seq/results/01_edgeRobject.Rds"))
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
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] DESeq2_1.22.2               SummarizedExperiment_1.12.0
##  [3] DelayedArray_0.8.0          BiocParallel_1.16.6        
##  [5] matrixStats_0.54.0          Biobase_2.42.0             
##  [7] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
##  [9] IRanges_2.16.0              S4Vectors_0.20.1           
## [11] BiocGenerics_0.28.0         RColorBrewer_1.1-2         
## [13] ggplot2_3.1.1               edgeR_3.24.3               
## [15] limma_3.38.3               
## 
## loaded via a namespace (and not attached):
##  [1] bit64_0.9-7            splines_3.5.1          Formula_1.2-3         
##  [4] assertthat_0.2.1       latticeExtra_0.6-28    blob_1.1.1            
##  [7] GenomeInfoDbData_1.2.0 yaml_2.2.0             RSQLite_2.1.1         
## [10] pillar_1.3.1           backports_1.1.4        lattice_0.20-35       
## [13] glue_1.3.1             digest_0.6.18          XVector_0.22.0        
## [16] checkmate_1.9.1        colorspace_1.4-1       htmltools_0.3.6       
## [19] Matrix_1.2-14          plyr_1.8.4             XML_3.98-1.15         
## [22] pkgconfig_2.0.2        genefilter_1.64.0      zlibbioc_1.28.0       
## [25] purrr_0.3.2            xtable_1.8-3           scales_1.0.0          
## [28] tibble_2.1.1           htmlTable_1.13.1       annotate_1.60.1       
## [31] withr_2.1.2            nnet_7.3-12            lazyeval_0.2.2        
## [34] survival_2.42-6        magrittr_1.5           crayon_1.3.4          
## [37] memoise_1.1.0          evaluate_0.13          foreign_0.8-71        
## [40] tools_3.5.1            data.table_1.12.2      stringr_1.4.0         
## [43] munsell_0.5.0          locfit_1.5-9.1         cluster_2.0.7-1       
## [46] AnnotationDbi_1.44.0   compiler_3.5.1         rlang_0.3.4           
## [49] RCurl_1.95-4.12        rstudioapi_0.10        htmlwidgets_1.3       
## [52] labeling_0.3           bitops_1.0-6           base64enc_0.1-3       
## [55] rmarkdown_1.12         gtable_0.3.0           DBI_1.0.0             
## [58] R6_2.4.0               gridExtra_2.3          knitr_1.22            
## [61] dplyr_0.8.0.1          bit_1.1-14             Hmisc_4.2-0           
## [64] stringi_1.4.3          Rcpp_1.0.1             geneplotter_1.60.0    
## [67] rpart_4.1-13           acepack_1.4.1          tidyselect_0.2.5      
## [70] xfun_0.6
```

