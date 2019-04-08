#####
## script to run CENTIPEDE for a given TF
## run RScript runCentipede.R TF
## where TF is the protein name
#####

#################
### FUNCTIONS ###
#################

createMotifMatches_inPeakSet <- function(TFname = "CTCF", peakSet){
  matches <- read.table(paste0(dir, "ATAC-seq/motifs/PWMScan_HOCOMOCOv10_mm10/", TFname, "_TFBS.bed"), stringsAsFactors = FALSE)
  matches <- GRanges(matches$V1, IRanges(start=matches$V2, end=matches$V3), 
                     strand=matches$V6, motif=matches$V4, score=matches$V5)
  ## keep those that overlap peaks
  matches <- matches[matches %over% peakSet]
  matches <- as.data.frame(matches)
  matches <- data.frame(motif_id=paste0(TFname, 1:nrow(matches)), 
                        'sequence-name'=matches$seqnames, start=matches$start, stop=matches$end, strand=matches$strand, 
                        score=matches$score, 'p-value'=0, 'q-value'=0, matched_sequence=matches$motif, stringsAsFactors = FALSE)
  colnames(matches)[2] <- "sequence-name"
  colnames(matches)[7] <- "p-value"
  colnames(matches)[8] <- "q-value"
  write.table(matches, file=paste0(dir, "ATAC-seq/motifs/", TFname, ".txt"), quote = FALSE, sep="\t", row.names = FALSE)
  print(paste0("Created motif binding sites overlapping peakSet for ", TFname))
  return(nchar(matches[1,9]))
}

insertionCountMatrix <- function(samples, TFname="CTCF"){
  ## use the function provided in the package CENTIPEDE.tutorial to create a matrix of insertion counts for each TF
  ## the expected input is the BAM file plus the output of FIMO with matches of the TF PWM along the genome
  ## we instead use 422 mouse TF with in-silico predicted TFBS based on the HOCOMOCO 10 database and PWMScan for mm10 (obtained from diffTF [https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz]; run with cutoff p-value 0.00001, background base composition 0.29;0.21;0.21;0.29)
  ## we retain all predicted binding sites, so we set the p-value filter to 1.
  data <- list()
  for(sample in samples){
    data[[sample]] <- centipede_data(bam_file = paste0(dir, "ATAC-seq/data/BWA/", sample, ".noDUPs.GQ.bam"),
                                     fimo_file = paste0(dir, "ATAC-seq/motifs/", TFname, ".txt"),
                                     pvalue = 1, flank_size = 100 )
  }
  return(data)
}

runCentipede <- function(insertionMatrix){
  fit <- list()
  for(sample in names(insertionMatrix)){
    fit[[sample]] <- fitCentipede( Xlist = list(ATACseq = insertionMatrix[[sample]]$mat), 
                                   Y = as.matrix(data.frame( Intercept = rep(1, nrow(insertionMatrix[[sample]]$mat)))),
                                   DampLambda=0.01,DampNegBin=0.001)
  }
  return(fit)
}



#################
library(csaw)
library(CENTIPEDE)
library(CENTIPEDE.tutorial)

dir <- "/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/"

args <- commandArgs(trailingOnly = TRUE)

## metadata
meta <- read.table(paste0(dir, "ATAC-seq/data/metadata_ATACseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
meta <- meta[meta$QCpass==1,]

## data
filtered.data <- readRDS(paste0(dir, "ATAC-seq/results/03_windowCounts_filteredWindows.Rds"))
merged <- mergeWindows(rowRanges(filtered.data), tol=150, max.width = 1500)
peakSet <- merged$region



# create binding sites file
motifLength <- createMotifMatches_inPeakSet(args[1], peakSet = peakSet)

# compute insertion-count matrix for all samples
counts <- insertionCountMatrix(samples = meta$sample, TFname = args[1])

# fit CENTIPEDE
fit <- runCentipede(insertionMatrix = counts)

# number of TFBSs per sample
probs <- mapply(`[[`, fit, 1)
row.names(probs) <- paste0(counts[[1]]$regions[,1], ":", counts[[1]]$regions[,2], "-", counts[[1]]$regions[,3])
write.table(probs, paste0(dir, "ATAC-seq/results/05_TFBSs/", args[1], "_probBindingSites.tab"), quote = FALSE, sep="\t")


# plot footprint
pdf(paste0(dir, "ATAC-seq/results/05_footprints/", args[1], ".pdf", width=10, height=10))
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(sample in names(fit)){
  plot(fit[[sample]]$LambdaParList$ATACseq[1:c(200+motifLength)], ylim=c(0,round(max(fit[[sample]]$LambdaParList$ATACseq[1:c(200+motifLength)]),3)), type="l", lwd=3, col="blue", 
       main=sample, ylab="Cut-site probability", xlab="distance to motif (bp)", axes=FALSE)
  lines(fit[[sample]]$LambdaParList$ATACseq[c(200+motifLength+1):c(400+motifLength*2)], lwd=3, col="red")
  box(bty="l"); axis(2)
  axis(1, at=c(0,50,100,120,170,220), labels=c(-100,-50,0,0,50,100))
  abline(v=c(100,120),lty=2)
  #legend("topright", legend = c("forward", "reverse"), col=c("blue", "red"), lwd=2)
}
dev.off()

saveRDS(fit, paste0(dir, "ATAC-seq/results/05_footprints/", args[1], "_centipedeFit.Rds"))


sessionInfo()

