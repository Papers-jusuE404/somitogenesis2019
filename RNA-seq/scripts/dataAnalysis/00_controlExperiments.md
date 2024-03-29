---
title: "Control experiment to ensure left and right somites are equivalent"
date: "07 August 2019"
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



Before creating matched RNA and ATAC-seq datasets using the two somites from each pair it was important to determine that the two somites from a pair are truly equivalent.

To test this, we collected the first two pairs of somites from two different embryos from the same litter. All somites were used to produce RNA-seq libraries. Additionally, the RNA from a somite pair was used to produce libraries again, but using only a quarter of the recommended reagents ('bis' samples); this was to assess whether using less reagents resulted in the same quality data.


```r
## metadata
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq_CONTROLlibraries.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta
```

```
##              sample library embryo somite side stage     date size
## 1      ctrl_e1_SI-1 do16433     e1     SI    1 20-25 03.08.17  434
## 2      ctrl_e1_SI-2 do16434     e1     SI    2 20-25 03.08.17  654
## 3     ctrl_e1_SII-1 do16435     e1    SII    1 20-25 03.08.17  509
## 4     ctrl_e1_SII-2 do16436     e1    SII    2 20-25 03.08.17  420
## 5      ctrl_e2_SI-1 do16437     e2     SI    1 20-25 03.08.17  432
## 6      ctrl_e2_SI-2 do16438     e2     SI    2 20-25 03.08.17  393
## 7     ctrl_e2_SII-1 do16439     e2    SII    1 20-25 03.08.17  461
## 8     ctrl_e2_SII-2 do16440     e2    SII    2 20-25 03.08.17  570
## 9  ctrl_e1_SII-1bis do16441     e1    SII    1 20-25 03.08.17  682
## 10 ctrl_e1_SII-2bis do16442     e1    SII    2 20-25 03.08.17  630
##        conc
## 1  12.51092
## 2  11.37400
## 3  12.28483
## 4  12.92770
## 5  17.52016
## 6  13.98273
## 7  12.27829
## 8  11.43267
## 9  34.50733
## 10 16.40458
```

### Quality control 

The ten samples were sequenced across one lane of an Illumina 4000, producing single-end 50bp reads. All samples were sequenced successfully, producing around 37+-7.1 million reads


```r
## counts
data <- read.table(paste0(dir, "RNA-seq/data/geneCounts_CONTROLlibraries.RAW.tsv"))
stopifnot(identical(meta$library, colnames(data)[-1]))

## mapping statistics (from STAR logs)
mapping.stats <- read.table(paste0(dir, "RNA-seq/data/mappingStatistics_CONTROLlibraries.tsv"), stringsAsFactors = FALSE, header = TRUE)
stopifnot(identical(colnames(data)[-1], mapping.stats$sample))

## counting statistics (from STAR gene counts)
counting.stats <- read.table(paste0(dir, "RNA-seq/data/countingStatistics_CONTROLlibraries.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(data)[-1], row.names(counting.stats)))
counting.stats$uniquelyMapped <- colSums(data[,-1])

ggplot(mapping.stats, aes(1, total/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("library size (million fragments)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/data-1.png)<!-- -->

From these, a median of 75.37% were uniquely mapped, plus 23.56% multimapped, which is expected for 50bp reads. Most of the uniquely mapped reads could be unambiguously assigned to annotated exons (85.44%).


```r
# summary(mapping.stats$unique/mapping.stats$total*100)
# summary(mapping.stats$multimapped/mapping.stats$total*100)

ggplot(counting.stats, aes(1, uniquelyMapped/1e6)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + ylab("total fragments uniquely mapped to exons (millions)") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/inExons-1.png)<!-- -->

```r
# summary(counting.stats$uniquelyMapped/mapping.stats$unique*100)
```

All samples but one had very uniform mapping statistics. The outlier (`do16440`) has higher proportion of multimapped reads (30.2%) which results in a lower proportion of uniquely mapped reads (67.5%). But from those, the proportion within exons is equivalent to other samples.

Similarly, all samples have a very uniform number of genes detected, around 22 thousand,


```r
counting.stats$nGenes <- apply(data[,-1], 2, function(x) sum(x>0))
# summary(counting.stats$nGenes)

ggplot(counting.stats, aes(1, nGenes/1e3)) + geom_violin(trim=FALSE) + geom_jitter(aes(1,nGenes/1e3, col=mapping.stats$total/1e6), width = 0.02) + ylab("total genes detected (thousands)") + scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "PiYG"))) + labs(colour="libsize") + xlab("") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![](00_controlExperiments_files/figure-html/numGenes-1.png)<!-- -->

Based on this, all samples seem to be of good quality.

### Normalisation

We normalise to account for differences in sequencing depth, using `edgeR`'s method.
We filter out lowly expressed genes:


```r
## create an edgeR object
meta$group <- paste(meta$embryo, meta$somite, sep=".")
y <- DGEList(counts=data[,-1], samples = meta, genes = data[,1], group = meta$group)

## filter low abundance genes
means <- aveLogCPM(y)
keep <- filterByExpr(y)

plot(density(means[keep]), lwd=2, xlab="average log counts-per-million", main="", bty="l", ylim=c(0,1), xlim=c(-5,15))
lines(density(means[-keep]), lty=2, lwd=2)
legend("topright", legend = c("kept", "filtered"), lty=c(1,2), lwd=2)
```

![](00_controlExperiments_files/figure-html/filter-1.png)<!-- -->

```r
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   38692   16844
```

```r
y <- y[keep, , keep.lib.sizes=FALSE]
```

And estimate normalisation factors, which successfully unify the samples.


```r
## normalisation
y <- calcNormFactors(y)

dataNorm <- cpm(y, log=TRUE, prior.count = 1)
## add gene names
dataNorm <- cbind(data[row.names(dataNorm),1], as.data.frame(dataNorm))
colnames(dataNorm)[1] <- "gene"
write.table(dataNorm, paste0(dir, "RNA-seq/data/geneCounts_CONTROLlibraries.NORM_logCPM.tsv"), quote = FALSE, sep="\t")

par(mfrow=c(1,2))
boxplot(log2(data[keep,-1]+1), main="RAW", ylab=expression('log'[2]*' counts + 1'), axes=FALSE); box(bty="l"); axis(2, las=2)
boxplot(dataNorm[,-1], main="NORMALISED", ylab=expression('log'[2]*' counts-per-million + 1'), axes=FALSE); box(bty="l"); axis(2, las=2)
```

![](00_controlExperiments_files/figure-html/norm-1.png)<!-- -->

A PCA on the normalised expression of the 1000 most variable genes reveals that the sample that was an outlier on mapping statistics also is an outlier in the PCA. Embryo of origin is the next largest source of variation.


```r
data.vst <- vst(as.matrix(data[keep,-1]))

vars <- rowVars(data.vst)
names(vars) <- row.names(data.vst)
pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)
# plot(prop.var, ylab="proportion of variance explained", xlab="PC", bty="l")

df <- as.data.frame(pca$x)
ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta$embryo), shape=meta$somite)) + labs(col="stage", shape="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
```

![](00_controlExperiments_files/figure-html/PCA-1.png)<!-- -->

Based on this, we remove this sample from downstream analyses.


```r
## remove aoutlier
outlier <- mapping.stats$sample[which.min(mapping.stats$unique/mapping.stats$total)]
meta <- meta[-which(meta$library==outlier),]
dataNorm <- dataNorm[,-which(colnames(dataNorm)==outlier)]

pca <- prcomp(t(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1]))

eigs <- pca$sdev^2
prop.var <- eigs / sum(eigs)

df <- as.data.frame(pca$x)
ggplot(df, aes(PC1, PC2)) + geom_point(aes(col=as.factor(meta$embryo), shape=meta$somite)) + labs(col="stage", shape="somite") + xlab(paste("PC1 - ", round(prop.var[1]*100,2), "% var explained")) + ylab(paste("PC2 - ", round(prop.var[2]*100,2), "% var explained"))
```

![](00_controlExperiments_files/figure-html/PCA_noOut-1.png)<!-- -->

### Differential expression analysis

We can use differential expression analysis to test whether the two somites from each pair have significant differences in gene expression. We use `edgeR` to estimate the dispersion.


```r
## design
design <- model.matrix(~0+somite+side, meta)

## remove outlier
y$samples <- y$samples[-which(y$samples$library == outlier),]
y$counts <- y$counts[,-which(colnames(y$counts)==outlier)]

## dispersion
y <- estimateDisp(y, design, robust = TRUE)
par(mfrow=c(1,2))
plotBCV(y)

## fit model
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

![](00_controlExperiments_files/figure-html/edger-1.png)<!-- -->

```r
# summary(fit$df.prior)
```

And we test for differences between side 1 and 2, while controlling for which somite was used. No significant genes are identified.


```r
## test
side <- glmQLFTest(fit, ncol(design))

## results
side <- as.data.frame(topTags(side, n=Inf))
# sum(side$FDR < 0.05)  # 0
write.table(side, paste0(dir, "RNA-seq/results/00_DEresults_side_allSomites.tsv"), quote = FALSE, sep="\t")

plot(side$logCPM, side$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side$FDR < 0.05, "red", "black"))
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/somiteTrios_perStage-1.png)<!-- -->

Alternatively, we can test separately the somiteI and somiteII pairs, using the two embryos as biological replicates (for somite II we only have one replicate for one of the sides). Again, no significant genes are identified.


```r
## somite I only
y.SI <- DGEList(counts=data[,c(2:3,6:7)], genes=data[,1], samples=meta[c(1:2,5:6),], group=meta$group[c(1:2,5:6)])
y.SI <- y.SI[keep, , keep.lib.sizes=FALSE]
y.SI <- calcNormFactors(y.SI)

design.SI <- model.matrix(~side, meta[c(1:2,5:6),])
y.SI <- estimateDisp(y.SI, design.SI, robust = TRUE)
fit.SI <- glmQLFit(y.SI, design.SI, robust=TRUE)
side.SI <- glmQLFTest(fit.SI, ncol(design.SI))
side.SI <- as.data.frame(topTags(side.SI, n=Inf))

par(mfrow=c(1,2))
plot(side.SI$logCPM, side.SI$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SI$FDR < 0.05, "red", "black"), main="somite I")
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))

## somite II only
y.SII <- DGEList(counts=data[,c(4:5,8)], genes=data[,1], samples=meta[c(3:4,7),], group=meta$group[c(3:4,7)])
y.SII <- y.SII[keep, , keep.lib.sizes=FALSE]
y.SII <- calcNormFactors(y.SII)

design.SII <- model.matrix(~side, meta[c(3:4,7),])
y.SII <- estimateDisp(y.SII, design.SII, robust = TRUE)
fit.SII <- glmQLFit(y.SII, design.SII, robust=TRUE)
side.SII <- glmQLFTest(fit.SII, ncol(design.SII))
side.SII <- as.data.frame(topTags(side.SII, n=Inf))

plot(side.SII$logCPM, side.SII$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SII$FDR < 0.05, "red", "black"), main="somite II")
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/perSomite-1.png)<!-- -->



Even though this experiment is not very well powered, any systematic clear differences between the two somites should be picked up. Thus, we can use the two somites as equivalent samples, to performed matched RNA and ATAC experiments.

### Differences between reagent concentration

For one somite pair, the library for sequencing was prepared twice, using a quarter of the recommended reagents. Since there are no significant differences between the left and right somites from the same pair, we can use them as replicates to compare the reagent concentrations. It's not ideal but good enough for a sanity check. If reducing the reagents led to systematic differences in a group of genes (e.g., those expressed at low levels, or with high/low GC content) we should see them.

We don't. No genes are significantly DE between the full and quartered reagents libraries.


```r
y.reagents <- DGEList(counts=data[,c(4:5,9:10)], genes=data[,1], samples=meta[c(3:4,8:9),], group=meta$group[c(3:4,8:9)])
y.reagents <- y.reagents[keep, , keep.lib.sizes=FALSE]
y.reagents <- calcNormFactors(y.reagents)

design.reagents <- model.matrix(~c(0,0,1,1))
y.reagents <- estimateDisp(y.reagents, design.reagents, robust = TRUE)
fit.reagents <- glmQLFit(y.reagents, design.reagents, robust=TRUE)
reagents <- glmQLFTest(fit.reagents, ncol(design.reagents))
reagents <- as.data.frame(topTags(reagents, n=Inf))

plot(reagents$logCPM, reagents$logFC, pch=16, cex=0.7, xlab="log2 mean expression", ylab="log2 fold-change", col=ifelse(side.SI$FDR < 0.05, "red", "black"), main="somite I")
abline(h=0, col="red", lwd=2)
legend("topright", legend = c("DE", "non-DE"), pch=16, col=c("red", "black"))
```

![](00_controlExperiments_files/figure-html/reagents-1.png)<!-- -->

And indeed the expression estimates are pretty well matched, as expected for technical replicates.


```r
par(mfrow=c(1,2))
plot(dataNorm[,'do16435'], dataNorm[,'do16441'], pch=16, cex=0.5, main="e1-SII-side1", xlab="full", ylab="quarter", bty="l")
abline(0,1,col="red", lwd=2)
plot(dataNorm[,'do16436'], dataNorm[,'do16442'], pch=16, cex=0.5, main="e1-SII-side2", xlab="full", ylab="quarter", bty="l")
abline(0,1,col="red", lwd=2)
```

![](00_controlExperiments_files/figure-html/scatter-1.png)<!-- -->

### Conclusions

These control experiments show that it is valid to consider the left and right somites from the same pair as equivalent samples, with the same transcriptome. Also, that using a lower amount of reagents in the library prep step doesn't affect the library produced. 

This is illustrated below, with the correlation coefficients between the transcriptomes of all the samples. The technical replicates -that only differ in the amount of reagents used- have the highest coefficients, very close to 1; and these are significantly higher than the correlations between any other two samples (wilcoxon rank sum test p-value = 0.003175). 

The samples coming from the same somite-pair, that only differ on whether they are from the left or right side are the next highest correlations. However, these correlation can be as *low* as those between adjacent somites. Overall, the correlation coefficients decrease as the number of variables that are different between samples increases, which is expected.


```r
corr <- cor(dataNorm[,-1], method="spearman")
# corr <- cor(dataNorm[names(vars[order(vars, decreasing = TRUE)])[1:1000],-1])

## group the different correlations
df <- data.frame(corr=unique(as.vector(corr))[-1])
## by the similarity of the two samples involved
df$group <- c("side","somite","somite","embryo","embryo-side","embryo-somite","somite","somite","somite","somite","embryo-side","embryo","embryo-somite","somite","somite","side","embryo-somite","embryo-somite","embryo","replicate","side","embryo-somite","embryo-somite","embryo","side","replicate","side","somite","embryo-somite","embryo-somite","somite","embryo-somite","embryo-somite","embryo","embryo-side","side")
## technical 'replicates' are the most similar
## followed by samples that only differ by 'side'
## then those that differ by the 'somite' taken (both sides considered equal)
## then those from a different 'embryo' but same somite and side
## versus 'embryo-side' where embryo and side are different
##lastly, 'embryo-somite' where both embryo and somite are different
df$group <- factor(as.character(df$group), levels=c("replicate", "side", "somite", "embryo", "embryo-side", "embryo-somite"))
ggplot(df, aes(group, corr)) + geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
```

![](00_controlExperiments_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
# wilcox.test(df[df$group=="replicate",]$corr, df[df$group!="replicate",]$corr) # 0.003175

# replicate corrs:  0.9856984 0.9805222 pearson
#                   0.9900903 0.9878620 spearman
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
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
##  [4] assertthat_0.2.1       statmod_1.4.30         latticeExtra_0.6-28   
##  [7] blob_1.1.1             GenomeInfoDbData_1.2.0 yaml_2.2.0            
## [10] RSQLite_2.1.1          pillar_1.3.1           backports_1.1.4       
## [13] lattice_0.20-35        glue_1.3.1             digest_0.6.18         
## [16] XVector_0.22.0         checkmate_1.9.1        colorspace_1.4-1      
## [19] htmltools_0.3.6        Matrix_1.2-14          plyr_1.8.4            
## [22] XML_3.98-1.15          pkgconfig_2.0.2        genefilter_1.64.0     
## [25] zlibbioc_1.28.0        purrr_0.3.2            xtable_1.8-3          
## [28] scales_1.0.0           tibble_2.1.1           htmlTable_1.13.1      
## [31] annotate_1.60.1        withr_2.1.2            nnet_7.3-12           
## [34] lazyeval_0.2.2         survival_2.42-6        magrittr_1.5          
## [37] crayon_1.3.4           memoise_1.1.0          evaluate_0.13         
## [40] foreign_0.8-71         tools_3.5.1            data.table_1.12.2     
## [43] stringr_1.4.0          munsell_0.5.0          locfit_1.5-9.1        
## [46] cluster_2.0.7-1        AnnotationDbi_1.44.0   compiler_3.5.1        
## [49] rlang_0.3.4            grid_3.5.1             RCurl_1.95-4.12       
## [52] rstudioapi_0.10        htmlwidgets_1.3        labeling_0.3          
## [55] bitops_1.0-6           base64enc_0.1-3        rmarkdown_1.12        
## [58] gtable_0.3.0           DBI_1.0.0              R6_2.4.0              
## [61] gridExtra_2.3          knitr_1.22             dplyr_0.8.0.1         
## [64] bit_1.1-14             Hmisc_4.2-0            stringi_1.4.3         
## [67] Rcpp_1.0.1             geneplotter_1.60.0     rpart_4.1-13          
## [70] acepack_1.4.1          tidyselect_0.2.5       xfun_0.6
```

