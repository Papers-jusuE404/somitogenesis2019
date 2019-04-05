---
title: "Quality control of somite ATAC-seq data"
date: "26 March 2019"
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



We have ATAC-seq data for 75 samples comprising somite trios at six different developmental stages. We have processed the data and retained only good-quality, unique alignments. 


```r
## metadata
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), stringsAsFactors = FALSE, header = TRUE)

## BAM files
bam.files <- paste0(dir, "ATACseq/data/BWA/", meta$sample, ".noDUPs.GQ.bam")
```

We now need to decide which samples are of enough good quality for downstream analyses.

We can use a number of metrics to assess this, including the insert size distribution (which should show the expected nucleosomal pattern), signal-to-noise ratio inferred from the ability to call peaks and enrichment of signal at the transcription start site (TSS) of genes.

## Insert size distribution

We have analysed the insert size distribution and scored the samples based on the strength of their nucleosomal pattern, from 0 to 4. Two thirds of the samples have a score of 2 or higher, representing samples with a decent to good nucleosomal pattern.


```r
colors <- c("firebrick3", "tomato2", "cornflowerblue", "royalblue3", "royalblue4")
names(colors) <- 0:4
pie(table(meta$insSizeDist), col=colors, main="insert-size distribution score")
```

![](02_qualityControl_files/figure-html/insSize-1.png)<!-- -->

Luckily, the samples with different scores are spread across stages, showing no systematic biases to a particular time-point. Stage 27 is biased to lower scores but all other stages are pretty balanced.


```r
barplot(t(table(meta$stage, meta$insSizeDist)/rowSums(table(meta$stage, meta$insSizeDist)))*100, col=colors, xlim=c(0,7), width = 0.75, space = 0.25, ylab="% samples in each score group")
legend("topright", legend = 0:4, col=colors, pch=15)
```

![](02_qualityControl_files/figure-html/scoreStage-1.png)<!-- -->

## Signal-to-noise ratio

Regions of open chromatin should show and enrichment of sequencing fragments and thus appear as *peaks*, or pileups of mapped fragments that are over background noise levels. One way to assess whether such enrichment is present in the data is to call peaks. We have used MACS2 for this, in each individual sample. 

Samples with good nucleosomal patterns tend to have higher number of peaks.


```r
blacklist <- readRDS(paste0(dir, "ATAC-seq/data/mm10-blacklist.rds"))
peaks <- list()

for(sample in meta$sample){
  peaks[[sample]] <- read.table(paste0(dir, "ATAC-seq/peaks/individualReps/", sample, "_peaks.broadPeak"))[,-6]
  colnames(peaks[[sample]]) <- c("chr", "start", "end", "name", "score", "fc", "pval", "qval")
  peaks[[sample]] <- GRanges(peaks[[sample]]$chr, IRanges(peaks[[sample]]$start, peaks[[sample]]$end), fc=peaks[[sample]]$fc, fdr=peaks[[sample]]$qval, score=peaks[[sample]]$score)
  remove <- unique(queryHits(findOverlaps(peaks[[sample]], blacklist)))
  peaks[[sample]] <- peaks[[sample]][-remove]
}

stopifnot(identical(meta$sample, names(peaks)))
meta$nPeaks <- unlist(lapply(peaks, length))

bpInPeaks <- unlist(lapply(peaks, function(x) sum(width(x))))

par(mfrow=c(1,2))
boxplot(meta$nPeaks/1e3~meta$insSizeDist, col=colors, xlab="nucleosomal quality score", ylab="number of peaks x 1000", axes=FALSE)
box(bty="l"); axis(1, at=1:5, labels=0:4); axis(2, las=2)
boxplot(bpInPeaks/1e6~meta$insSizeDist, col=colors, xlab="nucleosomal quality score", ylab="basepairs in peaks (millions)", axes=FALSE)
box(bty="l"); axis(1, at=1:5, labels=0:4); axis(2, las=2)
```

![](02_qualityControl_files/figure-html/peaks-1.png)<!-- -->

The majority of samples with a score of 0 have very low numbers of peaks, suggesting we have mostly captured background noise. Samples with scores of 2 or higher generally do well at peak calling, with 20 to 50 thousand total peaks. Samples with a score of 1 are in between, and contain a mixture of failed and decent to good peak calls. Nonetheless, all score groups have some samples that completely fail peak calling, despite having good nucleosomal patterns. Even without considering outliers, the interquartile range is quite large, with a two-fold difference in the total number of peaks between samples with good nucleosomal patterns.


```r
tmp <- sapply(0:4, function(x) summary(meta[meta$insSizeDist==x,]$nPeaks))
colnames(tmp) <- 0:4
tmp
```

```
##                 0        1        2        3        4
## Min.       41.000  5266.00   114.00    29.00  2792.00
## 1st Qu.  1558.250 15613.00 18709.75 20463.25 27343.50
## Median   5890.500 20978.50 36141.00 34885.50 51915.00
## Mean     9916.688 20637.00 33227.60 33489.46 40992.43
## 3rd Qu. 15632.750 27796.75 48914.00 50071.25 56209.00
## Max.    30981.000 35612.00 68145.00 60019.00 65135.00
```

A cutoff of 15 thousand peaks removes the lowest third of the data.

----

The *fraction of reads in peaks* is a good measurement of the signal-to-noise ratio. This will of course be correlated to the total number of peaks called and therefore to the quality of the nucleosomal pattern.

We use `bedtools` to count the number of fragments mapped across called peaks.


```r
## write the blacklist as a BED file to use with bedtools
write.table(as.data.frame(blacklist)[,1:3], file=paste0(dir, "ATAC-seq/data/mm10-blacklist.bed"), quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

frip <- rep(0, nrow(meta))
names(frip) <- meta$sample
for(sample in meta$sample){
  frip[sample] <- sum(read.table(paste0(dir, "ATAC-seq/peaks/individualReps/", sample, "_peaks.FRiP.bed"))[,10])/2
}

meta$readsInPeaks <- frip

boxplot(meta$readsInPeaks/meta$goodQuality*100~meta$insSizeDist, col=colors, ylab="FRiP", xlab="nucleosomal pattern score")
abline(h=3, lty=2)
```

![](02_qualityControl_files/figure-html/frip-1.png)<!-- -->

Samples with a score of 2 or higher have around 5-25% of reads in peaks. A cutoff of 3% removes the lowest third of the data.

### Reproducibility of peaks

Finally, if the called peaks represent true signal, we should find the same regions in many samples. To confirm this we check how many of the bases covered by one sample's peaks are also covered in the other samples' peak calls. Below, each sample is on the x-axis; 1 corresponds to the comparison with itself, and the boxplots indicate the distribution of proportion of overlap with the other 74 samples.


```r
## compute number of bp that overlap between peak sets
peakOverlap <- matrix(ncol=length(peaks), nrow=length(peaks))
for(i in 1:length(peaks)){
  for(j in 1:length(peaks)){
    peakOverlap[i,j] <- sum(width(suppressWarnings(intersect(peaks[[i]], peaks[[j]]))))
  }
}
row.names(peakOverlap) <- names(peaks)
colnames(peakOverlap) <- names(peaks)

## convert into proportion
tmp <- peakOverlap/diag(peakOverlap)  # diagonal contains overlap with itself (i.e. total bp in peaks)

o <- order(meta$nPeaks, decreasing = TRUE)
boxplot(tmp[,o], las=2, col=colors[meta$insSizeDist[o]+1], axes=FALSE, ylab="% of shared bp between peak sets")
box(bty="l"); axis(2, las=2); mtext(side=1, line=0.5, text = "samples by decreasing number of peaks")
```

![](02_qualityControl_files/figure-html/peakOverlap-1.png)<!-- -->

Overall, samples with scores of 2 or more tend to to have 70-95% of their bases shared with other samples. There are instead a group of samples, predominantly with score of 0 that share very little of their peaks with other samples, suggesting that those peaks are false positives.


## Transcritpion-start site enrichment

Another measure of the quality of an ATAC-seq library is the degree of enrichment of reads around genes' transcription start sites (TSSs). To check this, we use the RNA-seq data to retain the genes expressed at moderate to high levels (8,684 genes) and extract the coordinate of their most 5' TSS. We then write a BED file with the 2kb interval centred at the TSS.


```r
## retrieve TSS coordinates for all genes
## we used Ensembl annotation version 93 for the RNA-seq analysis. Use the same for consistentcy.
ensembl <- useMart(host = 'http://jul2018.archive.ensembl.org', 
                   biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

tss <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'transcription_start_site'), mart = ensembl)

## retain only expressed genes
expr <- read.table(paste0(dir, "RNA-seq/data/geneCounts_SOMITES.RAW.tsv"))[-c(1:4),]
expr <- expr[,-grep("PSM", colnames(expr))] ## we remove the PSM samples to use only the somites transcriptomes
# normalise
sf <- estimateSizeFactorsForMatrix(expr[,-1])
expr <- t(t(expr[,-1])/sf)
means <- rowMeans(expr)
keep <- means > 100 # use only moderatly to highly expressed genes
# summary(keep) # 8684  
genes <- names(means[keep])

## keep the earliest TSS for each expressed gene
tss <- tss[tss$ensembl_gene_id %in% genes,]
tmp <- tss[tss$ensembl_gene_id %in% unique(tss[duplicated(tss$ensembl_gene_id),1]),]
tss <- tss[!(tss$ensembl_gene_id %in% unique(tss[duplicated(tss$ensembl_gene_id),1])),]

min <- sapply(unique(tmp$ensembl_gene_id), function(x) min(tmp[tmp$ensembl_gene_id==x,4]))
tmp <- unique(tmp[,1:3])
tmp$transcription_start_site <- min[tmp$ensembl_gene_id]

tss <- rbind(tss, tmp)
row.names(tss) <- tss$ensembl_gene_id; tss$ensembl_gene_id <- NULL

## export as BED file (start is 0-based in BED files)
tmp <- data.frame(chr=factor(paste0("chr",tss$chromosome_name), levels=paste0("chr", c(1:19,'X'))), start=tss$transcription_start_site-1001, end=tss$transcription_start_site+1000, gene=tss$external_gene_name, id=row.names(tss))  # providing levels for chromosome allows us to sort in the right order; non-canonical chromosomes become NA
tmp <- tmp[!is.na(tmp$chr),]
tmp <- tmp[order(tmp$chr, tmp$start),]

write.table(tmp, file = paste0(dir, "ATAC-seq/data/TSScoords_expressedGenes_flanking1kb.bed"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
```

Using `betools` we can then count the number of Tn5 insertions at each basepair in these 2kb windows (see `ATAC-seq/scripts/dataProcessing/03.0_shiftingReads.md` for details). An insertion is regarded as the 5' and 3' ends of a paired-end fragment, after shifting the coordinates to account for the 9bp used by Tn5 for the insertion.

The aggregate profile for all genes gives us a measure of the number of insertions at each basepair, relative to the TSS. We use the first and last 100bp of the 2kb interval as a measurement of the background insertion level, and compute the fold-change of each bp against this. Thus, the flanks should start at 1 and the fold-change should increase if there is an excess of insertions as we approach the TSS.

We indeed see this, with values as high as 10-fold enrichment at the TSS for some samples.


```r
## read insertion counts
sample <- meta$sample[1]
tss <- list()
for(sample in meta$sample){
  tss[[sample]] <- matrix(read.table(paste0(dir, "ATAC-seq/TSS/", sample, ".TSScount.bed"), stringsAsFactors = FALSE)[,7], ncol = 2001, byrow = TRUE)
}
saveRDS(tss, paste0(dir, "ATAC-seq/results/02_TSSinsertionCounts.Rds"))

## normalise to background
tss.norm <- t(do.call("cbind", lapply(tss, function(x) colMeans(x)/mean(colMeans(x[,c(1:100,1901:2001)])))))

plot(density(tss.norm[,1001]), main="TSS enrichment across samples", bty="l")
abline(v=4, lty=2)
```

![](02_qualityControl_files/figure-html/countInsertionsRes-1.png)<!-- -->

Again, there is great variability in the degree of enrichment. We define a threshold to remove the tail of samples with low values; a cutoff of 4 discards the lowest third of the data.

Below are the enrichment profiles for each sample, ordered by decreasing number of peaks and coloured by nucleosomal profile score.


```r
par(mfrow=c(13,6), mar=c(2,2,2,2))
for(sample in meta[order(meta$nPeaks, decreasing = TRUE),]$sample){
  plot(rollmean(tss.norm[sample,], k=25), type="l", lwd=5, main=sample, xlab="", ylab="", col=colors[meta[meta$sample==sample,]$insSizeDist+1], ylim=c(0,10), axes=FALSE)
  box(bty="l"); axis(1, at=c(0,1000,2000), labels = c("-1kb","TSS","1kb")); axis(2,las=2)
  abline(h=4, lty=2, lwd=2)
}

stopifnot(identical(row.names(tss.norm),meta$sample))
meta$TSSscore <- tss.norm[,1001]

## save the metadata that has new columns
write.table(meta, paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
```

![](02_qualityControl_files/figure-html/plotTSSenrichment-1.png)<!-- -->

Again, there is a correlation of good TSS enrichment for samples with high numbers of peaks and good nucleosomal profile, while the bad-quality samples also have very poor TSS enrichment.

-----

As an aside, we observe the greatest enrichment at the TSS, but progressively smaller peaks at ~200bp intervals there on. This reflects the signal coming from nucleosome-bound genes. Above I've plotted the aggregate profile of all genes, regardless of strand. If the analysis is split for genes in each strand, mirror profiles occur, with the nucleosomal pattern downstream of the TSS. 

In the aggregate profiles we only observe the signal from genes in the + strand because signal for these genes is stronger (not sure why).


```r
ann <- read.table(paste0(dir, "RNA-seq/data/Mus_musculus.GRCm38.93.ann"), row.names = 1)
tmp <- read.table(paste0(dir, "ATAC-seq/data/TSScoords_expressedGenes_flanking1kb.bed"), stringsAsFactors = FALSE)
strands <- ann[tmp$V5,5]

par(mfrow=c(1,2))
plot(colSums(tss[['e17_SI-2']][strands=="+",]), ylab="insertion counts", xlab="", type="l", axes=FALSE, main="genes in + strand")
box(bty="l"); axis(1, at=c(0,1000,2000), labels = c("-1kb", "TSS", "1kb")); axis(2)
abline(v=seq(1000,1800,200), lty=2, col="grey")
plot(colSums(tss[['e17_SI-2']][strands=="-",]), ylab="insertion counts", xlab="", type="l", axes=FALSE, main="genes in - strand")
box(bty="l"); axis(1, at=c(0,1000,2000), labels = c("1kb", "TSS", "-1kb")); axis(2)
abline(v=seq(200,1000,200), lty=2, col="grey")
```

![](02_qualityControl_files/figure-html/strandedEnrichment-1.png)<!-- -->

## Quality control

Together, these different metrics allow us to identify the good samples that have great signal-to-noise ratio; the bad samples that are basically just background noise and need to be discarded; and the decent samples, that sit in between.

To assign a *pass* on each metric we require:

- Nucleosomal pattern score >= 2.

- Number of peaks > 15,000.

- FRiP > 3%.

- TSS score > 4.

And with these criteria, 67% of the samples have 3 or 4 passes.


```r
qc <- data.frame(nuclosome = ifelse(meta$insSizeDist>=2, 1, 0), nPeaks = ifelse(meta$nPeaks > 15000, 1, 0), frip = ifelse(meta$readsInPeaks/meta$goodQuality*100 >= 3, 1, 0), tss = ifelse(meta$TSSscore > 4, 1, 0))
row.names(qc) <- meta$sample
# round(table(rowSums(qc))/nrow(qc)*100,2)
# 0     1     2     3     4 
# 18.67 12.00  2.67 14.67 52.00 
upset(qc, mainbar.y.label = "number of samples", sets.x.label = "number of samples\nthat pass", text.scale=1.5)
```

![](02_qualityControl_files/figure-html/qc-1.png)<!-- -->

Thus, we retain 50 samples and get rid of the 25 samples that only pass two of the four criteria above. With this remaining set of samples we have generally 2-4 replicates per somite-stage combination, except for stage27 that looses more and SIII of stage35.


```r
stopifnot(identical(row.names(qc), meta$sample))
meta$QCpass <- ifelse(rowSums(qc)>2, 1, 0)
table(meta[meta$QCpass==1,]$stage, meta[meta$QCpass==1,]$somite)
```

```
##     
##      SI SII SIII
##   8   2   2    2
##   18  3   4    3
##   21  2   4    1
##   25  4   4    3
##   27  3   1    1
##   35  5   4    2
```

```r
write.table(meta, paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
```

Nonetheless, these should be enough samples for a decent analysis.



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
##  [1] biomaRt_2.34.2              ggplot2_3.0.0              
##  [3] UpSetR_1.3.3                zoo_1.8-3                  
##  [5] DESeq2_1.18.1               GenomicAlignments_1.18.1   
##  [7] SummarizedExperiment_1.12.0 DelayedArray_0.4.1         
##  [9] matrixStats_0.54.0          Biobase_2.38.0             
## [11] Rsamtools_1.34.1            Biostrings_2.50.2          
## [13] XVector_0.22.0              GenomicRanges_1.34.0       
## [15] GenomeInfoDb_1.18.2         IRanges_2.16.0             
## [17] S4Vectors_0.20.1            BiocGenerics_0.28.0        
## [19] RColorBrewer_1.1-2         
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.3.1             bit64_0.9-7            splines_3.5.1         
##  [4] Formula_1.2-3          assertthat_0.2.0       latticeExtra_0.6-28   
##  [7] blob_1.1.1             GenomeInfoDbData_1.0.0 progress_1.2.0        
## [10] yaml_2.2.0             RSQLite_2.1.0          pillar_1.3.0          
## [13] backports_1.1.2        lattice_0.20-35        glue_1.3.0            
## [16] digest_0.6.15          checkmate_1.9.1        colorspace_1.3-2      
## [19] htmltools_0.3.6        Matrix_1.2-14          plyr_1.8.4            
## [22] XML_3.98-1.15          pkgconfig_2.0.2        genefilter_1.60.0     
## [25] zlibbioc_1.24.0        purrr_0.2.5            xtable_1.8-2          
## [28] scales_1.0.0           BiocParallel_1.12.0    tibble_1.4.2          
## [31] htmlTable_1.12         annotate_1.56.2        withr_2.1.2           
## [34] nnet_7.3-12            lazyeval_0.2.1         survival_2.42-6       
## [37] magrittr_1.5           crayon_1.3.4           memoise_1.1.0         
## [40] evaluate_0.11          foreign_0.8-71         prettyunits_1.0.2     
## [43] tools_3.5.1            data.table_1.11.4      hms_0.4.2             
## [46] stringr_1.3.1          locfit_1.5-9.1         munsell_0.5.0         
## [49] cluster_2.0.7-1        AnnotationDbi_1.44.0   bindrcpp_0.2.2        
## [52] compiler_3.5.1         rlang_0.2.1            grid_3.5.1            
## [55] RCurl_1.95-4.11        rstudioapi_0.7         htmlwidgets_1.2       
## [58] labeling_0.3           bitops_1.0-6           base64enc_0.1-3       
## [61] rmarkdown_1.10         gtable_0.2.0           curl_3.2              
## [64] DBI_1.0.0              R6_2.2.2               gridExtra_2.3         
## [67] knitr_1.20             dplyr_0.7.6            bit_1.1-14            
## [70] bindr_0.1.1            Hmisc_4.1-1            rprojroot_1.3-2       
## [73] stringi_1.2.4          Rcpp_0.12.18           geneplotter_1.56.0    
## [76] rpart_4.1-13           acepack_1.4.1          tidyselect_0.2.4
```
