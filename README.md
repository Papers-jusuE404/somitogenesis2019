# Regulatory control of mouse somitogenesis

This repository contains the scripts used to process and analyse the genomic data for the mouse somitogenesis project.

## ATAC-seq data

The `scripts` directory contains:

- `dataProcessing` - all commands used to manipulate sequencing data files. It includes data pre-processing to clean the BAM files (remove PCR duplicates, low quality alignments, etc.); read shifting to obtain the insertion site coordinates; generation of BigWig files with insertion counts per basepair for all peaks, for visualisation.

- `dataAnalysis` - all R markdown documents with the code (and associated results) for normalisation, differential accessibility analysis, etc.

The `data` directory contains metadata text files and the BAM files used in analyses (these are not included in this repo).


## RNA-seq data

...