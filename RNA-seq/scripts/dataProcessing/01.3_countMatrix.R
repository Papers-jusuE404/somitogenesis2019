## generate count matrix from RNA-seq data
dir <- "/user01/group_folders/Personal/Ximena/SOMITES/somitogenesis2019/"

ann <- read.table(paste0(dir, "RNA-seq/data/Mus_musculus.GRCm38.96.ann"), row.names=1, stringsAsFactors = FALSE)
colnames(ann) <- c("gene", "chr", "start", "end", "strand", "biotype")

samples <- list.dirs(paste0(dir,"RNA-seq/data/STAR"), recursive=FALSE)

# with the first sample, add gene names
data <- read.table(paste0(samples[1], "/ReadsPerGene.out.tab"), row.names = 1, stringsAsFactors = FALSE)
data$gene <- ann[match(row.names(data), row.names(ann)),1]
data <- data[,c(4,1)]

# append counts for all others
for(i in 2:length(samples)){
  data <- cbind(data, read.table(paste0(samples[i], "/ReadsPerGene.out.tab"), row.names = 1, stringsAsFactors = FALSE)[,1])
}
colnames(data)[-1] <- unlist(lapply(samples, function(x) unlist(strsplit(x, "/"))[11]))

## edit column names to have sample names instead of library id
meta <- read.table(paste0(dir, "RNA-seq/data/metadata_RNAseq.tsv"), stringsAsFactors = FALSE, header = TRUE)
data <- data[,c("gene",intersect(colnames(data), meta$library))] ## remove samples from stage39 which are not replicated and will not be used
colnames(data)[-1] <- meta[match(colnames(data)[-1], meta$library),]$sample

## separate mapping stats
mapping.stats <- data[1:4,-1]
data <- data[-c(1:4),]

## save results
write.table(data, paste0(dir, "RNA-seq/data/geneCounts.RAW.tsv"), quote = FALSE, sep="\t")
write.table(t(mapping.stats), paste0(dir, "RNA-seq/data/mappingStatistics.tsv"), quote = FALSE, sep="\t")
