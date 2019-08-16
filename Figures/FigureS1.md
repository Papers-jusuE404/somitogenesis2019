---
title: "Figure S1"
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



### Figure S1

To show that the left and right somites from the same pair are equivalent. Differential expression analysis of four pairs of somites, comparing the left versus right somites.


```r
res <- read.table(paste0(dir, "RNA-seq/results/00_DEresults_side_allSomites.tsv"))

# pdf(paste0(dir, "Figures/PDFs/FigureS4_MAplot_sides.pdf"), width = 5, height = 5, useDingbats=FALSE)
plot(res$logCPM, res$logFC, pch=16, cex=0.7, bty="l", xlab=expression('log'[2]*' average expression'), ylab=expression('log'[2]*' fold-change'))
abline(h=0, lwd=2, col="red")
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](FigureS1_files/figure-html/MAplot-1.png)<!-- -->

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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     htmltools_0.3.6
##  [5] yaml_2.2.0      Rcpp_1.0.1      stringi_1.4.3   rmarkdown_1.12 
##  [9] knitr_1.22      stringr_1.4.0   xfun_0.6        digest_0.6.18  
## [13] evaluate_0.13
```

