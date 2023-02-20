##Filter out all the transcripts which are artifacts
setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/DTU/")
load("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/mergelib2.RData")
dim(mergelib)
sqanti <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/sqanti/filter_lib_ML_MLresult_classification.txt")

sqanti$transcript <- gsub("transcript:","",sqanti$isoform)
trans_list <- as.data.frame(sqanti$transcript)
enst <- as.data.frame(trans_list$`sqanti$transcript`[grepl("ENST",trans_list$`sqanti$transcript`)==TRUE])
colnames(enst) <- "trans"

sqanti_isoform <- sqanti[sqanti$filter_result=="Isoform",]

isoform_list <- as.data.frame(sqanti_isoform$transcript)
isoform_enst <- isoform_list[grepl("ENST", isoform_list$`sqanti_isoform$transcript`)==TRUE,]


merge <- append(enst$trans, isoform_list$`sqanti_isoform$transcript`)
uniq <- as.data.frame(unique(merge))

##Get list of isoforms
mergelib <- as.data.frame(mergelib)
mergelib$transcript <- rownames(mergelib)
mergelib2 <- mergelib[mergelib$transcript %in% uniq$`unique(merge)`,]

print("end1")
##use this as input for DTU

mergelib2 <- mergelib2[,-25551] ##check this omg
dim(mergelib2)
print("end2")

save(mergelib2, file="mergelib_0.5_ENST.Rdata")

#table(mergelib2)
