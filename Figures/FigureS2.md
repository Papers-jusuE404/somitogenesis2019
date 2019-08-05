---
title: "Figure S2"
date: "05 July 2019"
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



### Figure S2

Figure on normalisation of ATAC-seq data.


```r
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq_GQ.tsv"), stringsAsFactors = FALSE, header = TRUE)
```

#### Raw data

MA plots show trended biases.

Plots using 10kb bins.


```r
background <- readRDS(paste0(dir, "ATAC-seq/results/03_backgroundCounts_10kbBins.Rds"))
adj.counts <- cpm(asDGEList(background), log=TRUE)

# pdf(paste0(dir, "Figures/PDFs/FigureS2_MAplot_10kb.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
par(mfrow=c(1,2), mar=c(2,2,2,2))
for (i in c(11,40)) {
  mval <- adj.counts[,1]-adj.counts[,i]
  smoothScatter(x=(adj.counts[,1]+adj.counts[,i])/2, y=mval, xlab="A", ylab="M", main=paste("1 vs", i), ylim=c(-4,4))
  abline(h=0, col="red", lwd=2)
}
```

![](FigureS2_files/figure-html/MAplot-1.png)<!-- -->

```r
# dev.off()
```

#### Normalised data

Efficiency bias normalisation doesn't help.

Plots using high-abundance 150bp windows.


```r
filtered.data <- readRDS(file=paste0(dir, "ATAC-seq/results/03_windowCounts_filteredWindows.Rds"))

# average counts
abval <- aveLogCPM(asDGEList(filtered.data))
o <- order(abval)

# raw counts
adjc <- log2(assay(filtered.data)+0.5)

# normalised counts
filtered.data <- normFactors(filtered.data, se.out = TRUE)
re.adjc <- cpm(asDGEList(filtered.data), log=TRUE)

# pdf(paste0(dir, "Figures/PDFs/FigureS2_efficiencyNormalisation.pdf"), useDingbats=FALSE, width = 7, height = 7)
par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in c(11,40)){
  mval <- adjc[,1]-adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, ylab="M", xlab="Average logCPM", main=paste("Raw 1 vs",i), ylim=c(-5,5))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
for(i in c(11,40)){
  mval <- re.adjc[,1]-re.adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, ylab="M", xlab="Average logCPM", main=paste("Normalised 1 vs",i), ylim=c(-5,5))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
```

![](FigureS2_files/figure-html/efficiency-1.png)<!-- -->

```r
# dev.off()
```

Trended normalisation is necessary to remove the bias. 


```r
# normalised counts - loess
re.adjc <- adjc - assay(filtered.data, "offset")/log(2)

# pdf(paste0(dir, "Figures/PDFs/FigureS2_trendedNormalisation.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
par(mfrow=c(1,2), mar=c(2,2,2,2))
# for(i in c(11,40)){
#   mval <- adjc[,1]-adjc[,i]
#   fit <- loessFit(x=abval, y=mval)
#   smoothScatter(abval, mval, ylab="M", xlab="Average logCPM", main=paste("Raw 1 vs",i))
#   lines(abval[o], fit$fitted[o], col="red")
#   abline(h=0)
# }
for(i in c(11,40)){
  mval <- re.adjc[,1]-re.adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, ylab="M", xlab="Average logCPM", main=paste("Normalised 1 vs",i))
  lines(abval[o], fit$fitted[o], col="red")
  abline(h=0)
}
```

![](FigureS2_files/figure-html/trended-1.png)<!-- -->

```r
# dev.off()
```

#### PCA


```r
vars <- rowVars(as.matrix(re.adjc))
tmp <- re.adjc[order(vars, decreasing=TRUE)[1:5000],]
pca <- prcomp(t(tmp))
df <- cbind(pca$x, meta)

# pdf(paste0(dir, "Figures/PDFs/FigureS2_PCA_norm.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(colour=as.factor(df$stage))) + labs(colour="stage")
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(colour=df$readsInPeakSet/df$goodQuality*100)) + scale_color_gradientn(colors=rev(brewer.pal(n=10,"RdYlBu"))) + labs(colour="FRiP")
multiplot(plotlist = plots, cols=2)
```

![](FigureS2_files/figure-html/pca-1.png)<!-- -->

```r
# dev.off()
```


```r
# pdf(paste0(dir, "Figures/PDFs/FigureS2_PCA_FRiP.pdf"), useDingbats=FALSE, width = 3.5, height = 3.5)
plot(df$PC1, df$readsInPeakSet/df$goodQuality*100, bty="l", xlab="PC1", ylab="FRiP", pch=16)
abline(lm(df$readsInPeakSet/df$goodQuality*100~df$PC1))
mtext(side=3, line=-1, text = paste("r =", round(cor(df$PC1, df$readsInPeakSet/df$goodQuality*100),2), "(p-value = 5.276e-06)"), adj=1)
```

![](FigureS2_files/figure-html/corr-1.png)<!-- -->

```r
# dev.off()
```



```r
pcs <- read.table(paste0(dir, "ATAC-seq/results/03_pcs_residuals.tab"))

meta$group <- factor(paste(meta$stage, meta$somite, sep="."))
design <- model.matrix(~0+group, meta)
colnames(design) <- paste0("stage",levels(meta$group))

norm.counts.corr.pca <- removeBatchEffect(re.adjc, design=design, covariates = pcs[,1:18])

vars <- rowVars(norm.counts.corr.pca)
tmp <- norm.counts.corr.pca[order(vars, decreasing=TRUE)[1:5000],]
pca <- prcomp(t(tmp))
df <- cbind(pca$x, meta)

# pdf(paste0(dir, "Figures/PDFs/FigureS2_PCA_norm_corrected.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(colour=as.factor(df$stage))) + labs(colour="stage")
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(colour=df$readsInPeakSet/df$goodQuality*100)) + scale_color_gradientn(colors=rev(brewer.pal(n=10,"RdYlBu"))) + labs(colour="FRiP")
multiplot(plotlist = plots, cols=2)
```

![](FigureS2_files/figure-html/removePCs-1.png)<!-- -->

```r
# dev.off()
```


```r
# pdf(paste0(dir, "Figures/PDFs/FigureS2_PCA_FRiP_clean.pdf"), useDingbats=FALSE, width = 3.5, height = 3.5)
plot(df$PC1, df$readsInPeakSet/df$goodQuality*100, bty="l", xlab="PC1", ylab="FRiP", pch=16)
abline(lm(df$readsInPeakSet/df$goodQuality*100~df$PC1))
mtext(side=3, line=-1, text = paste("r =", round(cor(df$PC1, df$readsInPeakSet/df$goodQuality*100),2), "(p-value = 0.099)"), adj=1)
```

![](FigureS2_files/figure-html/corrClean-1.png)<!-- -->

```r
# dev.off()
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
##  [1] RColorBrewer_1.1-2          edgeR_3.24.3               
##  [3] limma_3.38.3                ggplot2_3.1.1              
##  [5] csaw_1.16.1                 SummarizedExperiment_1.12.0
##  [7] DelayedArray_0.8.0          BiocParallel_1.16.6        
##  [9] matrixStats_0.54.0          Biobase_2.42.0             
## [11] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
## [13] IRanges_2.16.0              S4Vectors_0.20.1           
## [15] BiocGenerics_0.28.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1               locfit_1.5-9.1          
##  [3] lattice_0.20-35          prettyunits_1.0.2       
##  [5] Rsamtools_1.34.1         Biostrings_2.50.2       
##  [7] assertthat_0.2.1         digest_0.6.18           
##  [9] R6_2.4.0                 plyr_1.8.4              
## [11] RSQLite_2.1.1            evaluate_0.13           
## [13] httr_1.4.0               pillar_1.3.1            
## [15] zlibbioc_1.28.0          rlang_0.3.4             
## [17] GenomicFeatures_1.34.8   progress_1.2.0          
## [19] lazyeval_0.2.2           blob_1.1.1              
## [21] Matrix_1.2-14            rmarkdown_1.12          
## [23] labeling_0.3             stringr_1.4.0           
## [25] RCurl_1.95-4.12          bit_1.1-14              
## [27] biomaRt_2.38.0           munsell_0.5.0           
## [29] compiler_3.5.1           rtracklayer_1.42.2      
## [31] xfun_0.6                 pkgconfig_2.0.2         
## [33] htmltools_0.3.6          tidyselect_0.2.5        
## [35] tibble_2.1.1             GenomeInfoDbData_1.2.0  
## [37] XML_3.98-1.15            withr_2.1.2             
## [39] crayon_1.3.4             dplyr_0.8.0.1           
## [41] GenomicAlignments_1.18.1 bitops_1.0-6            
## [43] gtable_0.3.0             DBI_1.0.0               
## [45] magrittr_1.5             scales_1.0.0            
## [47] KernSmooth_2.23-15       stringi_1.4.3           
## [49] XVector_0.22.0           tools_3.5.1             
## [51] bit64_0.9-7              glue_1.3.1              
## [53] purrr_0.3.2              hms_0.4.2               
## [55] yaml_2.2.0               AnnotationDbi_1.44.0    
## [57] colorspace_1.4-1         memoise_1.1.0           
## [59] knitr_1.22
```

