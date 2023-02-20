#To generate merged transcript count matrix filtered by barcodes remaining after filtering via short-read analysis
library("DTUrtle")
library("BiocParallel")
library("GenomicRanges") 
library("Gviz")
library("rtracklayer") 
library("stageR")
library("tximport")
library("DESeq2")
library("stringr")
library("tidyr")
library("data.table")   
library("Matrix")
library("plyr")
library("dplyr")
library(stringr)
library("SingleCellExperiment")
library(methods) 
biocpar <- BiocParallel::MulticoreParam(1)
#BiocManager::valid()
setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/")
##import data
mergelib <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/mergelib.tsv")

#import gtf Annotation to get transcript to gene mapping

all_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/all_gtf.tsv", header=FALSE)

#write.table(all_gtf, file="all_gtf.tsv", sep="\t", row.names = FALSE, col.names = FALSE)
all_gtf$transcript_id <- gsub("transcript_id transcript:", "",all_gtf$V9)
all_gtf$transcript_id <- gsub(";.*", "",all_gtf$transcript_id)
all_gtf$gene_id <- gsub(".*gene_id ","",all_gtf$V9)
all_gtf$gene_id <- gsub(";","",all_gtf$gene_id)



rownames(mergelib) <- mergelib$transcript_id
mergelib <- subset(mergelib,select=-c(transcript_id))
fun1 <- function(x) {
x[is.na(x)] <- 0
x
}

mergelib <- apply(mergelib, 2, fun1)

save(mergelib,file="mergelib2.RData")
write.table(mergelib,file="mergelib2.tsv", col.names=TRUE, sep="\t", row.names=FALSE)

dim(mergelib) 

#Make list of genes with more than one transcript
txInfo <- all_gtf[,10:11]
txInfo <- unique(txInfo)
txInfo_sp <- split(txInfo, txInfo$gene_id)
txInfo_sp2 <- txInfo_sp[lapply(txInfo_sp, nrow) >1] ##remove genes with one transcript
txInfo_sp3 <- as.data.frame(do.call(rbind,txInfo_sp2)) ##flatten
rownames(txInfo_sp3) <- txInfo_sp3$isoform_id ##is this supposed to be transcript_id??
colnames(txInfo_sp3)[1] <- "isoform_id"
txInfo <- txInfo_sp3
#===================================

#==============================================

##Filter count matrix with genes with only one isoform

mergelib <- mergelib[which(
  rownames(mergelib) %in% txInfo_sp3$isoform_id), ]

save(mergelib, file="mergelib3.Rdata")
