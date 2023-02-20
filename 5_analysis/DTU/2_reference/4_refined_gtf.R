#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/")
all.gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/all_gtf.tsv", header=FALSE)

split_m <- split(all.gtf, all.gtf$V9)


fun10 <- function(x) {
  if(sum(x$V3=="transcript")>1){
    nrow <- NROW(x)
    dnom <- sum(x$V3=="transcript")
    half <- nrow/dnom
    x <- x[c(1:half),]
    
  }
  x
}


split_m2 <- lapply(split_m,fun10)
split_m2_flat <-  as.data.frame(do.call(rbind, split_m2))

##Test
#split2 <- split(split_m2_flat, split_m2_flat$V9)
#check_2 <- as.data.frame(sapply(lapply(split2,function(x) x$V3=="transcript"), sum))

split_m2_flat$V9 <- gsub("transcript:", "\"transcript:",split_m2_flat$V9)
split_m2_flat$V9 <- gsub("; gene_id ", "; gene_id \"gene:",split_m2_flat$V9)
split_m2_flat$V9 <- gsub(";", "\";",split_m2_flat$V9)

write.table(split_m2_flat, file="refined_gtf.gtf",sep="\t",col.names = FALSE, row.names=FALSE, quote=FALSE)
