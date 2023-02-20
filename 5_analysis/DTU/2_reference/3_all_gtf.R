#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/")

#if(!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager") 
#}
#BiocManager::install("satuRn")
#BiocManager::install("AnnotationHub")
#BiocManager::install("DEXSeq")
library(satuRn)
library(SummarizedExperiment)
library(AnnotationHub)
library(ensembldb)
library(edgeR)
library(ggplot2)
library(DEXSeq)
library(stageR)


#import gtf Annotation to get transcript to gene mapping

libA_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libA.gtf")
libB_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libB.gtf")
libC_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libC.gtf")
libD_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libD.gtf")
libE_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libE.gtf")
libF_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libF.gtf")
libG_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libG.gtf")
libH_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/corrected_libH.gtf")

all_gtf <- rbind(libA_gtf, libB_gtf, libC_gtf, libD_gtf, libE_gtf, libF_gtf, libG_gtf, libH_gtf)
write.table(all_gtf, file="all_gtf.tsv", sep="\t", row.names = FALSE, col.names = FALSE)