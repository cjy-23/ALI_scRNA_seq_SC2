#setwd("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/dturtle_round_ceiling_strand/")

#if(!requireNamespace("remotes", quietly = TRUE)){
#  install.packages("remotes")
#}
#remotes::install_github("TobiTekath/DTUrtle")
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
#library(stringr)
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

##Round to nearest 10 
#list <- c("libA_gtf", "libB_gtf", "libC_gtf", "libD_gtf", "libE_gtf", "libF_gtf", "libG_gtf", "libH_gtf")
#list <- c(libA_gtf, libB_gtf, libC_gtf, libD_gtf, libE_gtf, libF_gtf, libG_gtf, libH_gtf)

libA_gtf$V4 <- round_any(libA_gtf$V4, 10, ceiling)
libB_gtf$V4 <- round_any(libB_gtf$V4, 10, ceiling)
libC_gtf$V4 <- round_any(libC_gtf$V4, 10, ceiling)
libD_gtf$V4 <- round_any(libD_gtf$V4, 10, ceiling)
libE_gtf$V4 <- round_any(libE_gtf$V4, 10, ceiling)
libF_gtf$V4 <- round_any(libF_gtf$V4, 10, ceiling)
libG_gtf$V4 <- round_any(libG_gtf$V4, 10, ceiling)
libH_gtf$V4 <- round_any(libH_gtf$V4, 10, ceiling)


libA_gtf$V5 <- round_any(libA_gtf$V5, 10, ceiling)
libB_gtf$V5 <- round_any(libB_gtf$V5, 10, ceiling)
libC_gtf$V5 <- round_any(libC_gtf$V5, 10, ceiling)
libD_gtf$V5 <- round_any(libD_gtf$V5, 10, ceiling)
libE_gtf$V5 <- round_any(libE_gtf$V5, 10, ceiling)
libF_gtf$V5 <- round_any(libF_gtf$V5, 10, ceiling)
libG_gtf$V5 <- round_any(libG_gtf$V5, 10, ceiling)
libH_gtf$V5 <- round_any(libH_gtf$V5, 10, ceiling)


libA_gtf$lib <- "libA"
libB_gtf$lib <- "libB"
libC_gtf$lib <- "libC"
libD_gtf$lib <- "libD"
libE_gtf$lib <- "libE"
libF_gtf$lib <- "libF"
libG_gtf$lib <- "libG"
libH_gtf$lib <- "libH"

#Split data

libA <- split(libA_gtf, libA_gtf$V9)
libB <- split(libB_gtf, libB_gtf$V9)
libC <- split(libC_gtf, libC_gtf$V9)
libD <- split(libD_gtf, libD_gtf$V9)
libE <- split(libE_gtf, libE_gtf$V9)
libF <- split(libF_gtf, libF_gtf$V9)
libG <- split(libG_gtf, libG_gtf$V9)
libH <- split(libH_gtf, libH_gtf$V9)


#RBIND ##569868

bind <- append(libA,libB)
bind <- append(bind, libC)
bind <- append(bind, libD)
bind <- append(bind, libE)
bind <- append(bind, libF)
bind <- append(bind, libG)
bind <- append(bind, libH)


##Get rid of "transcript line"

bind2 <- lapply(bind, function(x) x[-1,])

##Find number of exons per transcript annotation
rows <- lapply(bind2, function(x) nrow(x))

##Find max number of exons ##71!!
rows <- as.data.frame(rows)
rows <- t(rows)
rows <- as.data.frame(rows)
max(rows$V1)

##make new nested data frame and then paste everything per nested dataframe (test run)

#try <- bind2$`transcript_id transcript:ENSG00000000003_100629025_100636732_3; gene_id ENSG00000000003;`
#chrom <- try[1,1]
#rest1 <- try[,4:5]
#rest1$tot <- paste(rest1$V4,rest1$V5,sep="_")
#rows <- nrow(rest1)
#rest2 <- paste(chrom,rest1$tot[1], rest1$tot[2],rest1$tot[3],rest1$tot[4],rest1$tot[5],rest1$tot[6],rest1$tot[7],rest1$tot[8], sep="_")
#cols <- c(1:rows)
#rest2 <- do.call(paste, c(rest1$tot, sep="_"))
#for (co in cols) data[co] <- NULL
#rest3 <- paste(chrom, rest1$tot[1:8],sep="_")

#gene <- sub(".*gene_id ", "", libA$`transcript_id transcript:ENSG00000000003_100629025_100636732_3; gene_id ENSG00000000003;`[1,9])

##make new nested data frame and then paste everything per nested dataframe (real)
df <- data.frame()
fun <- function(x) {
  chrom <- x[1,1]
  strand <- x[1,7]
  gene <- sub(".*gene_id ", "", x[1,9])
  gene <- sub(";","",gene)
  rest1 <- x[,4:5]
  rest1$tot <- paste(rest1$V4,rest1$V5,sep="_")
  #row <- nrow(rest1)
  rest2 <- paste(chrom,strand,gene,
                 rest1$tot[1], rest1$tot[2],rest1$tot[3],rest1$tot[4],rest1$tot[5],rest1$tot[6],rest1$tot[7],rest1$tot[8], rest1$tot[9], rest1$tot[10],
                 rest1$tot[11], rest1$tot[12],rest1$tot[13],rest1$tot[14],rest1$tot[15],rest1$tot[16],rest1$tot[17],rest1$tot[18], rest1$tot[19], rest1$tot[20],
                 rest1$tot[21], rest1$tot[22],rest1$tot[23],rest1$tot[24],rest1$tot[25],rest1$tot[26],rest1$tot[27],rest1$tot[28], rest1$tot[29], rest1$tot[30],
                 rest1$tot[31], rest1$tot[32],rest1$tot[33],rest1$tot[34],rest1$tot[35],rest1$tot[36],rest1$tot[37],rest1$tot[38], rest1$tot[39], rest1$tot[40],
                 rest1$tot[41], rest1$tot[42],rest1$tot[43],rest1$tot[44],rest1$tot[45],rest1$tot[46],rest1$tot[47],rest1$tot[48], rest1$tot[49], rest1$tot[50],
                 rest1$tot[51], rest1$tot[52],rest1$tot[53],rest1$tot[54],rest1$tot[55],rest1$tot[56],rest1$tot[57],rest1$tot[58], rest1$tot[59], rest1$tot[60],
                 rest1$tot[61], rest1$tot[62],rest1$tot[63],rest1$tot[64],rest1$tot[65],rest1$tot[66],rest1$tot[67],rest1$tot[68], rest1$tot[69], rest1$tot[70],
                 rest1$tot[71],sep="_")
  df <- rbind(df,rest2)
}

df <- lapply(bind2,fun)


##Flatten results

df_flat <- as.data.frame(unlist(df))

#nm <- as.data.frame(names(df))



##Duplicates
dup <- as.data.frame(duplicated(df_flat$`unlist(df)`))
dup_pos <- as.data.frame(which(duplicated(df_flat$`unlist(df)`))) ##which positions are duplicated? 
colnames(dup_pos) <- "name"
##=================add lib data

##make new nested data frame and then paste everything per nested dataframe (real)
df2 <- data.frame() ##569868
fun <- function(x) {
  chrom <- x[1,1]
  strand <- x[1,7]
  rest1 <- x[,4:5]
  gene <- sub(".*gene_id ", "", x[1,9])
  gene <- sub(";","",gene)
  rest1$tot <- paste(rest1$V4,rest1$V5,sep="_")
  lib <- x[1,10]
  #row <- nrow(rest1)
  rest2 <- paste(chrom,strand,gene,
                 rest1$tot[1], rest1$tot[2],rest1$tot[3],rest1$tot[4],rest1$tot[5],rest1$tot[6],rest1$tot[7],rest1$tot[8], rest1$tot[9], rest1$tot[10],
                 rest1$tot[11], rest1$tot[12],rest1$tot[13],rest1$tot[14],rest1$tot[15],rest1$tot[16],rest1$tot[17],rest1$tot[18], rest1$tot[19], rest1$tot[20],
                 rest1$tot[21], rest1$tot[22],rest1$tot[23],rest1$tot[24],rest1$tot[25],rest1$tot[26],rest1$tot[27],rest1$tot[28], rest1$tot[29], rest1$tot[30],
                 rest1$tot[31], rest1$tot[32],rest1$tot[33],rest1$tot[34],rest1$tot[35],rest1$tot[36],rest1$tot[37],rest1$tot[38], rest1$tot[39], rest1$tot[40],
                 rest1$tot[41], rest1$tot[42],rest1$tot[43],rest1$tot[44],rest1$tot[45],rest1$tot[46],rest1$tot[47],rest1$tot[48], rest1$tot[49], rest1$tot[50],
                 rest1$tot[51], rest1$tot[52],rest1$tot[53],rest1$tot[54],rest1$tot[55],rest1$tot[56],rest1$tot[57],rest1$tot[58], rest1$tot[59], rest1$tot[60],
                 rest1$tot[61], rest1$tot[62],rest1$tot[63],rest1$tot[64],rest1$tot[65],rest1$tot[66],rest1$tot[67],rest1$tot[68], rest1$tot[69], rest1$tot[70],
                 rest1$tot[71],lib, sep="_")
  df2 <- rbind(df2,rest2)
}


df2 <- lapply(bind2,fun)
nm <- names(df2)

##Flatten results

df2_flat <- as.data.frame(unlist(df2))
df2_flat$ID <- nm

##check good i
#====================duplicated data # 338,470
duplist <- as.data.frame(df2_flat[dup_pos$name,])

duplist1 <- as.data.frame(df2_flat[dup_pos$name,])

##switch oder
#duplist <- duplist[,c(3,1,2)]
## Split lib data

duplist[c('Transcript', 'Lib')] <- str_split_fixed(duplist$`unlist(df2)`, 'lib', 2)

  


##check good |
## Split data base =====================================

split <- split(duplist,duplist$Transcript)


##Split data all

df2_flat_2 <- df2_flat
df2_flat_2[c('Transcript', 'Lib')] <- str_split_fixed(df2_flat_2$`unlist(df2)`, 'lib', 2)

split_all <- split(df2_flat_2,df2_flat_2$Transcript)
split_all_backup <- split_all
#split_all <- split_all_backup


##Check matching (should be 231,398, actually)

#not_dup <- as.data.frame(df2_flat[!(df2_flat %in% duplist1),])
#not_dup <- as.data.frame(subset(df2_flat,!(df2_flat %in% duplist1)))
#not_dup <- setdiff(df2_flat, duplist1)
#not_dup1 <- as.data.frame(not_dup)



##====================make transcript dict and merge duplicates
##order based on ID decreasing so that ENST is preference
split_all <- lapply(split_all,function(x) x[order(x$ID, decreasing = TRUE),])

#tteststst <- split_all$ "chr3_150544290_150544500_150546050_150546500_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_NA_"

fun2 <- function(x) {
  x$after <- x[1,2]
  x
}

split_all_2 <- lapply(split_all,fun2)

##Check good |
##Flatten

flatlist <-  as.data.frame(do.call(rbind, split_all_2))
#colnames(flatlist)[2] <- "Before" 
#write.table(flatlist, file="flatlist_dict.tsv", row.names=TRUE,col.names = TRUE, sep="\t")


##moment of truth - check if there are any duplicate names in entries without duplicates

##subset non duplicate data

#nono <- as.data.frame(no_dup)
#nono <- NROW(df2_flat)
#nono_list <- 1:nono
#dup_pos2 <- as.list(dup_pos$name)
#diff <- setdiff(nono_list, dup_pos2)

#sub <- df2_flat[diff,]
#sub[c('Transcript', 'Lib')] <- str_split_fixed(sub$`unlist(df2)`, 'lib', 2)
#sub_check <- as.data.frame(sub$ID)
#sub_checking <- which(duplicated(sub_check))
#sub_checking_2 <- unique(sub_check)

##Check good |
#flatlist_split1 <- split(flatlist,flatlist$Transcript)
##================check no dups
fun3 <- function(x) {
  if(NROW(x)>1) {
  x<-NULL
  }  else {
    return(x)
  }
}

split_all_nodup <- lapply(split_all_2, fun3)

what <- lapply(split_all_nodup,function(x) is.null(x))
which <- as.data.frame(which(what==FALSE))

split_all_nodup2 <- split_all_2[which$`which(what == FALSE)`]


ID <- lapply(split_all_nodup2, function(x) return(x$after))
flat_ID <- as.data.frame(do.call(rbind,ID))

##check dup
sub_checking <- as.data.frame(which(duplicated(flat_ID$V1)))
sub_checking_2 <- unique(flat_ID)

duplicate <- as.data.frame(flat_ID[sub_checking$`which(duplicated(flat_ID$V1))`,])
colnames(duplicate) <- "dup"

flatlist2 <- flatlist[flatlist$after %in% duplicate$dup,]

##=======================================================================
##======================== take split and only take the first row (final names)
round1 <- lapply(split_all_2, function(x) x[1,])
round1_flat <- as.data.frame(do.call(rbind,round1))

sub_checkingg <- as.data.frame(which(duplicated(round1_flat$ID)))
#duplicatee <- as.data.frame(round1_flat[sub_checkingg$`which(duplicated(round1_flat$ID))`,])

#split - where are the duplicates

splitt <- split(round1_flat,round1_flat$after) ##Use ID instead of after because essentially the same

splitt_sub <- splitt[lapply(splitt,nrow)>1] ##Which ID's are duplicated but have different reference exon boundaries


#sort by order of lib
splitt_sub1 <- splitt_sub
#splitt_sub1 <- lapply(splitt_sub,function(x) x[order(x$Lib),])

##extract transcript ID 1
fun4 <- function(x){
  x$ID1 <- sub("\\;.*", "", x$ID)
  x
}
splitt_sub2 <- lapply(splitt_sub1,fun4)



##extract transcript ID 2
fun5 <- function(x){
  x$ID2 <- strsplit(x$ID1, "transcript_id transcript:")[[1]][2]
  x
}
splitt_sub3 <- lapply(splitt_sub2,fun5)


##check good |

##new column
fun6 <- function(x){
  x$new <- paste(x$ID2,seq(1,by = 1,length.out = nrow(x)),sep="_")
  x
}

splitt_sub4 <- lapply(splitt_sub3, fun6)


##Check good |

##flatten this new dict file

splitt_sub4_flatlist <-  as.data.frame(do.call(rbind, splitt_sub4))



##extract the after and new column

#int_dict <- splitt_sub4_flatlist[,c(7,8)]
#colnames(int_dict)[1] <- "ID2"

##redesign the flatlist according to int dict




##merge flatlist and splitt_sub4_flatlist

flatlist_3 <- flatlist

flatlist_3$ID1 <- sub("\\;.*", "", flatlist_3$after)
flatlist_3$ID2 <- sub(".*transcript_id transcript:", "", flatlist_3$ID1)

##Index matching
#mmm <- as.data.frame(merge(flatlist_3$ID2, int_dict[, c("ID2", "new")], by="ID2"))
                     
#DDD <- flatlist_3$ID2
#EEE <- int_dict$ID2
#FFF <- as.data.frame(int_dict$ID2 %in% flatlist_3$ID2)

#int_dict
#as.data.frame(which(!is.na(match(flatlist_3$ID2, int_dict$after))))

#flatlist_3$new <- ifelse(is.na(match(flatlist_3$ID2,int_dict$after))=="FALSE",int_dict$new,flatlist_3$ID2)

#merge(table1, table2[, c("pid", "val2")], by="pid")
#flatlist_3$new <- lapply(flatlist_3,fun7)


#flatlist_3_ID1 <- apply(flatlist_3_test,1,fun7)
#flatlist_3_ID1 <-as.data.frame(flatlist_3_ID1)


#flatlist_3_test <- flatlist[4]

#flatlist_3$ID1 <- lapply(flatlist_3, strsplit(flatlist_3$after, ";")[[1]][1])
#flatlist_3$ID2 <- strsplit(flatlist_3$ID1, "transcript_id transcript:")[[1]][2]

##make dict file

merge <- merge(flatlist_3, splitt_sub4_flatlist, by=c("ID2","Transcript"),all.x=TRUE)
merge_split <-split(merge,merge$Transcript)
#merge_split2 <- merge_split[lapply(merge_split,nrow)>1] 

##Don't really need this but just in case
fun9 <- function(x) {
    x$new <- x[1,13]
    x
  }

##check equal
merge_split2 <- lapply(merge_split,fun9) 
##check equal
identical(merge_split, merge_split2)

##flatten

merge_flat <-  as.data.frame(do.call(rbind, merge_split2))

merge_flat$new <- ifelse(is.na(merge_flat$new)=="TRUE",  merge_flat$ID2, merge_flat$new)

merge_flat$before1<- sub("\\;.*", "", merge_flat$ID.x)
merge_flat$before2 <- sub(".*transcript_id transcript:", "", merge_flat$before1)

final_dict <- merge_flat[,c(15,1,13,5)]


colnames(final_dict) <- c("Before", "Mid", "After", "Lib")
write.table(final_dict, file="final_dictionary_all.tsv", sep="\t", col.names=TRUE, row.names=TRUE)


#Split per library

lib_split <- split(final_dict, final_dict$Lib)

libA_dict <- lib_split$A
libB_dict <- lib_split$B
libC_dict <- lib_split$C
libD_dict <- lib_split$D
libE_dict <- lib_split$E
libF_dict <- lib_split$F
libG_dict <- lib_split$G
libH_dict <- lib_split$H

##Save file

write.table(libA_dict, file="libA_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libB_dict, file="libB_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libC_dict, file="libC_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libD_dict, file="libD_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libE_dict, file="libE_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libF_dict, file="libF_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libG_dict, file="libG_dict", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(libH_dict, file="libH_dict", sep="\t", col.names=TRUE, row.names=TRUE)


