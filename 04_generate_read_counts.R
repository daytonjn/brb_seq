#generate read/count matrix from starsolo

setwd("/cluster/tufts/dopmanlab/Jacob/brbseq/")
meta <- read.table("/cluster/tufts/dopmanlab/Jacob/brbseq/brbseq_combined.csv", header=TRUE,
                   sep = ',')
meta <- meta[,2:ncol(meta)]


library(data.table)
library(Matrix)

#reorder matrix counts for map_d3_v1
matrix_dir <- "./map_d3_v2_norrna/Solo.out/Gene/raw/"
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")

#f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
counts_d3 <- as.data.frame(as.matrix(readMM(f)))
close(f)

feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE,
                         stringsAsFactors = FALSE, data.table = F)
colnames(counts_d3) <- barcode.names$V1
rownames(counts_d3) <- feature.names$V1
rownames(meta) <- meta$code

#verify that all row labels in meta are columns in data
all(colnames(counts_d3) == meta$barcode)
idx = match(meta$barcode, colnames(counts_d3))
counts_d3 = counts_d3[,idx]
all(colnames(counts_d3) == meta$barcode)
colnames(counts_d3) <- meta$code
all(colnames(counts_d3) == meta$code)

fwrite(counts_d3, file = "umi.counts_d3_v2.txt", sep = "\t", quote = F, row.names = T,
       col.names = T)
###
#reorder matrix counts for map_d2_v1
matrix_dir <- "./map_d2_v2_norrna/Solo.out/Gene/raw/"

f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
#f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
counts_d2 <- as.data.frame(as.matrix(readMM(f)))
close(f)

feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE,
                      stringsAsFactors = FALSE, data.table = F)
colnames(counts_d2) <- barcode.names$V1
rownames(counts_d2) <- feature.names$V1
rownames(meta) <- meta$code

#verify that all row labels in meta are columns in data
all(colnames(counts_d2) == meta$barcode)
idx = match(meta$barcode, colnames(counts_d2))
counts_d2 = counts_d2[,idx]
all(colnames(counts_d2) == meta$barcode)
colnames(counts_d2) <- meta$code
all(colnames(counts_d2) == meta$code)

fwrite(counts_d2, file = "umi.counts_d2_v2.txt", sep = "\t", quote = F, row.names = T,
       col.names = T)
