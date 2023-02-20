#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/")

#if(!requireNamespace("remotes", quietly = TRUE)){
#  install.packages("remotes")
#}
#remotes::install_github("TobiTekath/round_ceiling_strand_gene")
#library("dturtle")
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
biocpar <- BiocParallel::MulticoreParam(1)


##import gtf

libA_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libA_12.4.gtf", header=FALSE)
libB_gtf  <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libB_12.4.gtf", header=FALSE)
libC_gtf  <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libC_12.4.gtf", header=FALSE)
libD_gtf  <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libD_12.4.gtf", header=FALSE)
libE_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libE_12.4.gtf", header=FALSE)
libF_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libF_12.4.gtf", header=FALSE)
libG_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libG_12.4.gtf", header=FALSE)
libH_gtf  <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/novel_isoforms/isoform_annotated.filtered_libH_12.4.gtf", header=FALSE)

#import lib_dicts

libA_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libA_dict")
libB_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libB_dict")
libC_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libC_dict")
libD_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libD_dict")
libE_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libE_dict")
libF_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libF_dict")
libG_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libG_dict")
libH_dict <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/libH_dict")

##==================================
##==================================

#Let's make new updated gtf files!
libA_gtf2 <- libA_gtf
libA_gtf2$before1<- sub("\\;.*", "", libA_gtf2$V9)
libA_gtf2$before2 <- sub(".*transcript_id transcript:", "", libA_gtf2$before1)
colnames(libA_gtf2)[11] <- "Before"
libA_gtf2$lol1 <- sub(".*gene_id", "", libA_gtf2$V9)
libA_gtf2$lol2 <- paste("gene_id",libA_gtf2$lol1, sep="")
##merge libA_dict

libA_gtf2_m <- merge(libA_gtf2, libA_dict, by="Before", all.x=TRUE)
libA_gtf2_m$new_V9 <- paste("transcript_id transcript:",libA_gtf2_m$After,";"," ",libA_gtf2_m$lol2, sep="")


##remove duplicates

libA_gt2_m_split <- split(libA_gtf2_m, libA_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libA_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libA_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libA_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libA_gt2_m_split_2 <- libA_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libA_gt2_m_split3 <- libA_gt2_m_split  
libA_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libA_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libA_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libA_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libA_gt2_m_split4 <- lapply(libA_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libA_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libA_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libA_gt2_m_split4))


libA_gt2_m_split4_flat$V9 <- libA_gt2_m_split4_flat$new_V9
rownames(libA_gt2_m_split4_flat) <- NULL
libA_gt2_m_split4_flat <- libA_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libA_gt2_m_split4_flat, file="corrected_libA.gtf",sep="\t",col.names = TRUE)

##libB
#Let's make new updated gtf files!
libB_gtf2 <- libB_gtf
libB_gtf2$before1<- sub("\\;.*", "", libB_gtf2$V9)
libB_gtf2$before2 <- sub(".*transcript_id transcript:", "", libB_gtf2$before1)
colnames(libB_gtf2)[11] <- "Before"
libB_gtf2$lol1 <- sub(".*gene_id", "", libB_gtf2$V9)
libB_gtf2$lol2 <- paste("gene_id",libB_gtf2$lol1, sep="")
##merge libB_dict

libB_gtf2_m <- merge(libB_gtf2, libB_dict, by="Before", all.x=TRUE)
libB_gtf2_m$new_V9 <- paste("transcript_id transcript:",libB_gtf2_m$After,";"," ",libB_gtf2_m$lol2, sep="")


##remove duplicates

libB_gt2_m_split <- split(libB_gtf2_m, libB_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libB_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libB_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libB_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libB_gt2_m_split_2 <- libB_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libB_gt2_m_split3 <- libB_gt2_m_split  
libB_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libB_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libB_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libB_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libB_gt2_m_split4 <- lapply(libB_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libB_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libB_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libB_gt2_m_split4))


libB_gt2_m_split4_flat$V9 <- libB_gt2_m_split4_flat$new_V9
rownames(libB_gt2_m_split4_flat) <- NULL
libB_gt2_m_split4_flat <- libB_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libB_gt2_m_split4_flat, file="corrected_libB.gtf",sep="\t",col.names = TRUE)



##libC
#Let's make new updated gtf files!
libC_gtf2 <- libC_gtf
libC_gtf2$before1<- sub("\\;.*", "", libC_gtf2$V9)
libC_gtf2$before2 <- sub(".*transcript_id transcript:", "", libC_gtf2$before1)
colnames(libC_gtf2)[11] <- "Before"
libC_gtf2$lol1 <- sub(".*gene_id", "", libC_gtf2$V9)
libC_gtf2$lol2 <- paste("gene_id",libC_gtf2$lol1, sep="")
##merge libC_dict

libC_gtf2_m <- merge(libC_gtf2, libC_dict, by="Before", all.x=TRUE)
libC_gtf2_m$new_V9 <- paste("transcript_id transcript:",libC_gtf2_m$After,";"," ",libC_gtf2_m$lol2, sep="")


##remove duplicates

libC_gt2_m_split <- split(libC_gtf2_m, libC_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libC_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libC_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libC_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libC_gt2_m_split_2 <- libC_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libC_gt2_m_split3 <- libC_gt2_m_split  
libC_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libC_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libC_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libC_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libC_gt2_m_split4 <- lapply(libC_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libC_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libC_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libC_gt2_m_split4))


libC_gt2_m_split4_flat$V9 <- libC_gt2_m_split4_flat$new_V9
rownames(libC_gt2_m_split4_flat) <- NULL
libC_gt2_m_split4_flat <- libC_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libC_gt2_m_split4_flat, file="corrected_libC.gtf",sep="\t",col.names = TRUE)

##libD
#Let's make new updated gtf files!
libD_gtf2 <- libD_gtf
libD_gtf2$before1<- sub("\\;.*", "", libD_gtf2$V9)
libD_gtf2$before2 <- sub(".*transcript_id transcript:", "", libD_gtf2$before1)
colnames(libD_gtf2)[11] <- "Before"
libD_gtf2$lol1 <- sub(".*gene_id", "", libD_gtf2$V9)
libD_gtf2$lol2 <- paste("gene_id",libD_gtf2$lol1, sep="")
##merge libD_dict

libD_gtf2_m <- merge(libD_gtf2, libD_dict, by="Before", all.x=TRUE)
libD_gtf2_m$new_V9 <- paste("transcript_id transcript:",libD_gtf2_m$After,";"," ",libD_gtf2_m$lol2, sep="")


##remove duplicates

libD_gt2_m_split <- split(libD_gtf2_m, libD_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libD_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libD_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libD_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libD_gt2_m_split_2 <- libD_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libD_gt2_m_split3 <- libD_gt2_m_split  
libD_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libD_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libD_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libD_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libD_gt2_m_split4 <- lapply(libD_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libD_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libD_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libD_gt2_m_split4))


libD_gt2_m_split4_flat$V9 <- libD_gt2_m_split4_flat$new_V9
rownames(libD_gt2_m_split4_flat) <- NULL
libD_gt2_m_split4_flat <- libD_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libD_gt2_m_split4_flat, file="corrected_libD.gtf",sep="\t",col.names = TRUE)

##libE
#Let's make new updated gtf files!
libE_gtf2 <- libE_gtf
libE_gtf2$before1<- sub("\\;.*", "", libE_gtf2$V9)
libE_gtf2$before2 <- sub(".*transcript_id transcript:", "", libE_gtf2$before1)
colnames(libE_gtf2)[11] <- "Before"
libE_gtf2$lol1 <- sub(".*gene_id", "", libE_gtf2$V9)
libE_gtf2$lol2 <- paste("gene_id",libE_gtf2$lol1, sep="")
##merge libE_dict

libE_gtf2_m <- merge(libE_gtf2, libE_dict, by="Before", all.x=TRUE)
libE_gtf2_m$new_V9 <- paste("transcript_id transcript:",libE_gtf2_m$After,";"," ",libE_gtf2_m$lol2, sep="")


##remove duplicates

libE_gt2_m_split <- split(libE_gtf2_m, libE_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libE_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libE_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libE_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libE_gt2_m_split_2 <- libE_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libE_gt2_m_split3 <- libE_gt2_m_split  
libE_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libE_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libE_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libE_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libE_gt2_m_split4 <- lapply(libE_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libE_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libE_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libE_gt2_m_split4))


libE_gt2_m_split4_flat$V9 <- libE_gt2_m_split4_flat$new_V9
rownames(libE_gt2_m_split4_flat) <- NULL
libE_gt2_m_split4_flat <- libE_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libE_gt2_m_split4_flat, file="corrected_libE.gtf",sep="\t",col.names = TRUE)

##libF
#Let's make new updated gtf files!
libF_gtf2 <- libF_gtf
libF_gtf2$before1<- sub("\\;.*", "", libF_gtf2$V9)
libF_gtf2$before2 <- sub(".*transcript_id transcript:", "", libF_gtf2$before1)
colnames(libF_gtf2)[11] <- "Before"
libF_gtf2$lol1 <- sub(".*gene_id", "", libF_gtf2$V9)
libF_gtf2$lol2 <- paste("gene_id",libF_gtf2$lol1, sep="")
##merge libF_dict

libF_gtf2_m <- merge(libF_gtf2, libF_dict, by="Before", all.x=TRUE)
libF_gtf2_m$new_V9 <- paste("transcript_id transcript:",libF_gtf2_m$After,";"," ",libF_gtf2_m$lol2, sep="")


##remove duplicates

libF_gt2_m_split <- split(libF_gtf2_m, libF_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libF_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libF_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libF_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libF_gt2_m_split_2 <- libF_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libF_gt2_m_split3 <- libF_gt2_m_split  
libF_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libF_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libF_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libF_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libF_gt2_m_split4 <- lapply(libF_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libF_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libF_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libF_gt2_m_split4))


libF_gt2_m_split4_flat$V9 <- libF_gt2_m_split4_flat$new_V9
rownames(libF_gt2_m_split4_flat) <- NULL
libF_gt2_m_split4_flat <- libF_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libF_gt2_m_split4_flat, file="corrected_libF.gtf",sep="\t",col.names = TRUE)

##libG
#Let's make new updated gtf files!
libG_gtf2 <- libG_gtf
libG_gtf2$before1<- sub("\\;.*", "", libG_gtf2$V9)
libG_gtf2$before2 <- sub(".*transcript_id transcript:", "", libG_gtf2$before1)
colnames(libG_gtf2)[11] <- "Before"
libG_gtf2$lol1 <- sub(".*gene_id", "", libG_gtf2$V9)
libG_gtf2$lol2 <- paste("gene_id",libG_gtf2$lol1, sep="")
##merge libG_dict

libG_gtf2_m <- merge(libG_gtf2, libG_dict, by="Before", all.x=TRUE)
libG_gtf2_m$new_V9 <- paste("transcript_id transcript:",libG_gtf2_m$After,";"," ",libG_gtf2_m$lol2, sep="")


##remove duplicates

libG_gt2_m_split <- split(libG_gtf2_m, libG_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libG_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libG_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libG_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libG_gt2_m_split_2 <- libG_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libG_gt2_m_split3 <- libG_gt2_m_split  
libG_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libG_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libG_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libG_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libG_gt2_m_split4 <- lapply(libG_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libG_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libG_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libG_gt2_m_split4))


libG_gt2_m_split4_flat$V9 <- libG_gt2_m_split4_flat$new_V9
rownames(libG_gt2_m_split4_flat) <- NULL
libG_gt2_m_split4_flat <- libG_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libG_gt2_m_split4_flat, file="corrected_libG.gtf",sep="\t",col.names = TRUE)


##libH
#Let's make new updated gtf files!
libH_gtf2 <- libH_gtf
libH_gtf2$before1<- sub("\\;.*", "", libH_gtf2$V9)
libH_gtf2$before2 <- sub(".*transcript_id transcript:", "", libH_gtf2$before1)
colnames(libH_gtf2)[11] <- "Before"
libH_gtf2$lol1 <- sub(".*gene_id", "", libH_gtf2$V9)
libH_gtf2$lol2 <- paste("gene_id",libH_gtf2$lol1, sep="")
##merge libH_dict

libH_gtf2_m <- merge(libH_gtf2, libH_dict, by="Before", all.x=TRUE)
libH_gtf2_m$new_V9 <- paste("transcript_id transcript:",libH_gtf2_m$After,";"," ",libH_gtf2_m$lol2, sep="")


##remove duplicates

libH_gt2_m_split <- split(libH_gtf2_m, libH_gtf2_m$new_V9)

check_2 <- as.data.frame(sapply(lapply(libH_gt2_m_split,function(x) x$V3=="transcript"), sum))
check_3 <- as.data.frame(check_2$`sapply(lapply(libH_gt2_m_split, function(x) x$V3 == "transcript"), sum)`>1)
colnames(check_3) <- "check"
##Check which position in libH_gt2_m_split


check_4 <- as.data.frame(which(check_3$check=="TRUE"))
libH_gt2_m_split_2 <- libH_gt2_m_split[check_4$`which(check_3$check == "TRUE")`]
check_5 <- as.data.frame(check_2[check_4$`which(check_3$check == "TRUE")`,])


##flatten
libH_gt2_m_split3 <- libH_gt2_m_split  
libH_gt2_m_split_2_flat <-  as.data.frame(do.call(rbind, libH_gt2_m_split_2))
check_6 <- as.data.frame(sapply(lapply(libH_gt2_m_split_2,function(x) x$V3=="exon"), sum))
check_7 <- as.data.frame(libH_gt2_m_split_2_flat$V3)



fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


libH_gt2_m_split4 <- lapply(libH_gt2_m_split3,fun10)


check_8 <- as.data.frame(sapply(lapply(libH_gt2_m_split4,function(x) x$V3=="transcript"), sum))

##let's flatten!!

libH_gt2_m_split4_flat <-  as.data.frame(do.call(rbind, libH_gt2_m_split4))


libH_gt2_m_split4_flat$V9 <- libH_gt2_m_split4_flat$new_V9
rownames(libH_gt2_m_split4_flat) <- NULL
libH_gt2_m_split4_flat <- libH_gt2_m_split4_flat[,-c(1,11:17)]



write.table(libH_gt2_m_split4_flat, file="corrected_libH.gtf",sep="\t",col.names = TRUE)














  
