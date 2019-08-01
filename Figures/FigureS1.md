---
title: "Figure S1"
date: "02 July 2019"
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



### Figure S1

This figure exemplifies the criteria used to QC the ATAC-seq data.


```r
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
```

#### Insert size distribution


```r
## in 01_insertSizeDistribution.Rmd we obtain the fragmen sizes for each clean BAM
diagnostics <- readRDS(paste0(dir, "ATAC-seq/results/01_mappingStats_ATACseq.Rds"))
stopifnot(identical(substr(names(diagnostics), 1, nchar(names(diagnostics))-14), meta$sample))

## get a representative sample from each insert size score group
samples <- c("e23_SII-2", "e16_SI-2", "e26_SII-2", "e32_SII-2", "e14_SII-2")

# pdf(paste0(dir, "Figures/PDFs/FigureS1_insSizeDist.pdf"), useDingbats=FALSE, width = 15, height = 3)
par(mfrow=c(1,5), mar=c(2,2,2,2))
for(i in 0:4){
  plot(density(diagnostics[[grep(samples[i+1], names(diagnostics))]]$sizes), axes=FALSE, main=paste("score =",i), xlab="", ylab="", lwd=3, xlim=c(0,1000))
  box(bty="l"); axis(1)
}
```

![](FigureS1_files/figure-html/insSize-1.png)<!-- -->

```r
# dev.off()
```


```r
plots <- list()
plots[[1]] <- ggplot(meta, aes(as.factor(insSizeDist), size)) + geom_boxplot() + ggtitle("DNA fragment size") + xlab("insert size dist score") + ylab("bp") + th
plots[[2]] <- ggplot(meta, aes(as.factor(insSizeDist), log10(conc))) + geom_boxplot() + ggtitle("DNA concentration") + xlab("insert size dist score") + ylab(expression('log'[10]*' [nM]')) + th
plots[[3]] <- ggplot(meta, aes(as.factor(insSizeDist), goodQuality/1e6)) + geom_boxplot() + ggtitle("Total unique GQ alignments") + xlab("insert size dist score") + ylab("million") + th

# pdf(paste0(dir, "Figures/PDFs/FigureS1_insSizeDist_corrs.pdf"), useDingbats=FALSE, width = 12, height = 4)
multiplot(plotlist = plots, cols=3)
```

![](FigureS1_files/figure-html/corrs-1.png)<!-- -->

```r
# dev.off()
```

#### TSS enrichment


```r
tss <- readRDS(paste0(dir, "ATAC-seq/results/02_TSSinsertionCounts.Rds"))
tss.norm <- t(do.call("cbind", lapply(tss, function(x) colMeans(x)/mean(colMeans(x[,c(1:100,1901:2001)])))))

# pdf(paste0(dir, "Figures/PDFs/FigureS1_TSSenrich.pdf"), useDingbats=FALSE, width = 15, height = 3)
par(mfrow=c(1,5), mar=c(2,2,2,2))
for(sample in samples){
  plot(rollmean(tss.norm[sample,], k=25), type="l", lwd=3, main=sample, xlab="", ylab="", ylim=c(1,8), axes=FALSE)
  box(bty="l"); axis(1, at=c(0,1000,2000), labels = c("-1kb","TSS","1kb")); axis(2,las=2)
  abline(h=5, lty=2, lwd=2)
}
```

![](FigureS1_files/figure-html/tss-1.png)<!-- -->

```r
# dev.off()
```

#### Peak calling


```r
# data_summary <- function(x) {
#    m <- mean(x)
#    ymin <- m-sd(x)
#    ymax <- m+sd(x)
#    return(c(y=m,ymin=ymin,ymax=ymax))
# }

## highlight the samples used above
meta$sel <- ifelse(meta$sample %in% samples, 1, 0)

plots <- list()
plots[[1]] <- ggplot(meta, aes(as.factor(insSizeDist), nPeaks/1e3, colour=sel)) + geom_violin() + geom_jitter(width=0.15) + xlab("insert size dist score") + ylab("number of peaks x 1000") + geom_hline(yintercept = 15, lty=2, col="grey") + th + theme(legend.position = "none")
plots[[2]] <- ggplot(meta, aes(as.factor(insSizeDist), readsInPeaks/goodQuality*100, colour=sel)) + geom_violin() + geom_jitter(width=0.15) +  xlab("insert size dist score") + ylab("FRiP") + geom_hline(yintercept = 3, lty=2, col="grey") + th + theme(legend.position = "none")

# pdf(paste0(dir, "Figures/PDFs/FigureS1_numberPeaks.pdf"), useDingbats=FALSE, width = 8, height = 4)
multiplot(plotlist = plots, cols=2)
```

![](FigureS1_files/figure-html/peakCalls-1.png)<!-- -->

```r
# dev.off()
```


#### Quality control


```r
qc <- data.frame(nuclosome = ifelse(meta$insSizeDist>=2, 1, 0), nPeaks = ifelse(meta$nPeaks > 15000, 1, 0), frip = ifelse(meta$readsInPeaks/meta$goodQuality*100 >= 3, 1, 0), tss = ifelse(meta$TSSscore > 4, 1, 0))
row.names(qc) <- meta$sample

# pdf(paste0(dir, "Figures/PDFs/FigureS1_upsetPlot.pdf"), useDingbats=FALSE, width = 5, height = 5)
upset(qc, mainbar.y.label = "number of samples", sets.x.label = "number of samples\nthat pass", text.scale=1.25)
```

![](FigureS1_files/figure-html/qc-1.png)<!-- -->

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
##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] UpSetR_1.3.3                zoo_1.8-5                  
##  [3] ggplot2_3.1.1               csaw_1.16.1                
##  [5] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
##  [7] BiocParallel_1.16.6         matrixStats_0.54.0         
##  [9] Biobase_2.42.0              Rsamtools_1.34.1           
## [11] Biostrings_2.50.2           XVector_0.22.0             
## [13] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
## [15] IRanges_2.16.0              S4Vectors_0.20.1           
## [17] BiocGenerics_0.28.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1               locfit_1.5-9.1          
##  [3] lattice_0.20-35          prettyunits_1.0.2       
##  [5] assertthat_0.2.1         digest_0.6.18           
##  [7] R6_2.4.0                 plyr_1.8.4              
##  [9] RSQLite_2.1.1            evaluate_0.13           
## [11] httr_1.4.0               pillar_1.3.1            
## [13] zlibbioc_1.28.0          rlang_0.3.4             
## [15] GenomicFeatures_1.34.8   progress_1.2.0          
## [17] lazyeval_0.2.2           blob_1.1.1              
## [19] Matrix_1.2-14            rmarkdown_1.12          
## [21] labeling_0.3             stringr_1.4.0           
## [23] RCurl_1.95-4.12          bit_1.1-14              
## [25] biomaRt_2.38.0           munsell_0.5.0           
## [27] compiler_3.5.1           rtracklayer_1.42.2      
## [29] xfun_0.6                 pkgconfig_2.0.2         
## [31] htmltools_0.3.6          tidyselect_0.2.5        
## [33] gridExtra_2.3            tibble_2.1.1            
## [35] GenomeInfoDbData_1.2.0   edgeR_3.24.3            
## [37] XML_3.98-1.15            withr_2.1.2             
## [39] crayon_1.3.4             dplyr_0.8.0.1           
## [41] GenomicAlignments_1.18.1 bitops_1.0-6            
## [43] gtable_0.3.0             DBI_1.0.0               
## [45] magrittr_1.5             scales_1.0.0            
## [47] stringi_1.4.3            limma_3.38.3            
## [49] tools_3.5.1              bit64_0.9-7             
## [51] glue_1.3.1               purrr_0.3.2             
## [53] hms_0.4.2                yaml_2.2.0              
## [55] AnnotationDbi_1.44.0     colorspace_1.4-1        
## [57] memoise_1.1.0            knitr_1.22
```

