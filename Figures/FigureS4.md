---
title: "Figure S4"
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



### Figure S4

Figure on batch effect removal from RNA-seq data.


```r
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$use==1,]

# data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"))
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"))
dataNorm <- read.table(paste0(dir, "RNA-seq/data/geneCounts.NORM_logCPM.tsv"))

pcs <- read.table(paste0(dir, "RNA-seq/results/02_pcs_residuals.tab"))
```

PCA of the normalised data, without controlling for batch effects. There is clear separation by batch within each stage.


```r
## variance-stabilise data
data <- data[row.names(dataNorm), colnames(dataNorm)]
data.vst <- vst(as.matrix(data[,-1]))

vars <- rowVars(data.vst)
names(vars) <- row.names(dataNorm)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$stage))) + labs(col="stage") + xlab("PC1") + ylab("PC2")
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$date))) + labs(col="batch") + xlab("PC1") + ylab("PC2")

# pdf(paste0(dir, "Figures/PDFs/FigureS3_noCorrection.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
multiplot(plotlist = plots, cols=2)
```

![](FigureS4_files/figure-html/norm-1.png)<!-- -->

```r
# dev.off()
```

After regressing out the first 14 PCs on the residuals of the fit of interest the batch effect is gone.


```r
dataNorm.pca <- removeBatchEffect(data.vst, design = model.matrix(~0+group, meta), covariates = pcs[,1:14])

vars <- rowVars(dataNorm.pca)
names(vars) <- row.names(dataNorm.pca)
pca <- prcomp(t(dataNorm.pca[names(vars[order(vars, decreasing = TRUE)])[1:1000],]))

df <- as.data.frame(pca$x)
plots <- list()
plots[[1]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$stage))) + labs(col="stage") + xlab("PC1") + ylab("PC2")
plots[[2]] <- ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta[meta$use==1,]$date))) + labs(col="batch") + xlab("PC1") + ylab("PC2")

# pdf(paste0(dir, "Figures/PDFs/FigureS3_correction.pdf"), useDingbats=FALSE, width = 7, height = 3.5)
multiplot(plotlist = plots, cols=2)
```

![](FigureS4_files/figure-html/norm.corrected-1.png)<!-- -->

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
##  [1] RColorBrewer_1.1-2          ggplot2_3.1.1              
##  [3] sva_3.30.1                  genefilter_1.64.0          
##  [5] mgcv_1.8-24                 nlme_3.1-137               
##  [7] DESeq2_1.22.2               SummarizedExperiment_1.12.0
##  [9] DelayedArray_0.8.0          BiocParallel_1.16.6        
## [11] matrixStats_0.54.0          Biobase_2.42.0             
## [13] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
## [15] IRanges_2.16.0              S4Vectors_0.20.1           
## [17] BiocGenerics_0.28.0         edgeR_3.24.3               
## [19] limma_3.38.3               
## 
## loaded via a namespace (and not attached):
##  [1] bit64_0.9-7            splines_3.5.1          Formula_1.2-3         
##  [4] assertthat_0.2.1       latticeExtra_0.6-28    blob_1.1.1            
##  [7] GenomeInfoDbData_1.2.0 yaml_2.2.0             pillar_1.3.1          
## [10] RSQLite_2.1.1          backports_1.1.4        lattice_0.20-35       
## [13] glue_1.3.1             digest_0.6.18          XVector_0.22.0        
## [16] checkmate_1.9.1        colorspace_1.4-1       htmltools_0.3.6       
## [19] Matrix_1.2-14          plyr_1.8.4             XML_3.98-1.15         
## [22] pkgconfig_2.0.2        zlibbioc_1.28.0        purrr_0.3.2           
## [25] xtable_1.8-3           scales_1.0.0           htmlTable_1.13.1      
## [28] tibble_2.1.1           annotate_1.60.1        withr_2.1.2           
## [31] nnet_7.3-12            lazyeval_0.2.2         survival_2.42-6       
## [34] magrittr_1.5           crayon_1.3.4           memoise_1.1.0         
## [37] evaluate_0.13          foreign_0.8-71         tools_3.5.1           
## [40] data.table_1.12.2      stringr_1.4.0          munsell_0.5.0         
## [43] locfit_1.5-9.1         cluster_2.0.7-1        AnnotationDbi_1.44.0  
## [46] compiler_3.5.1         rlang_0.3.4            RCurl_1.95-4.12       
## [49] rstudioapi_0.10        htmlwidgets_1.3        labeling_0.3          
## [52] bitops_1.0-6           base64enc_0.1-3        rmarkdown_1.12        
## [55] gtable_0.3.0           DBI_1.0.0              R6_2.4.0              
## [58] gridExtra_2.3          knitr_1.22             dplyr_0.8.0.1         
## [61] bit_1.1-14             Hmisc_4.2-0            stringi_1.4.3         
## [64] Rcpp_1.0.1             geneplotter_1.60.0     rpart_4.1-13          
## [67] acepack_1.4.1          tidyselect_0.2.5       xfun_0.6
```

