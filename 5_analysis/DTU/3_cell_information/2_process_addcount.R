library(Seurat)
library(dplyr)
library(purrr)
library(tibble)
library(dplyr)
library(Libra)
library(xlsx)
library(stringr)
setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/")

#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/")
libA_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libA_scp_meta.txt")
libA_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libA_match_cell_barcode")

libB_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libB_scp_meta.txt")
libB_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libB_match_cell_barcode")


libC_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libC_scp_meta.txt")
libC_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libC_match_cell_barcode")

libD_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libD_scp_meta.txt")
libD_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libD_match_cell_barcode")

libE_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libE_scp_meta.txt")
libE_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libE_match_cell_barcode")

libF_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libF_scp_meta.txt")
libF_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libF_match_cell_barcode")

libG_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libG_scp_meta.txt")
libG_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libG_match_cell_barcode")

libH_scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/libH_scp_meta.txt")
libH_match_cell_barcode <- read.csv("/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/flames_barcodes/libH_match_cell_barcode")



#-------------Add column in scp_meta.txt for count

colnames(libA_match_cell_barcode) <-c("NAME", "Count")
libA_scp_meta_2 <- merge(libA_scp_meta, libA_match_cell_barcode, by="NAME")

colnames(libB_match_cell_barcode) <-c("NAME", "Count")
libB_scp_meta_2 <- merge(libB_scp_meta, libB_match_cell_barcode, by="NAME")

colnames(libC_match_cell_barcode) <-c("NAME", "Count")
libC_scp_meta_2 <- merge(libC_scp_meta, libC_match_cell_barcode, by="NAME")


colnames(libD_match_cell_barcode) <-c("NAME", "Count")
libD_scp_meta_2 <- merge(libD_scp_meta, libD_match_cell_barcode, by="NAME")

colnames(libE_match_cell_barcode) <-c("NAME", "Count")
libE_scp_meta_2 <- merge(libE_scp_meta, libE_match_cell_barcode, by="NAME")

colnames(libF_match_cell_barcode) <-c("NAME", "Count")
libF_scp_meta_2 <- merge(libF_scp_meta, libF_match_cell_barcode, by="NAME")

colnames(libG_match_cell_barcode) <-c("NAME", "Count")
libG_scp_meta_2 <- merge(libG_scp_meta, libG_match_cell_barcode, by="NAME")


colnames(libH_match_cell_barcode) <-c("NAME", "Count")
libH_scp_meta_2 <- merge(libH_scp_meta, libH_match_cell_barcode, by="NAME")



#Save data
write.table(libA_scp_meta_2, file="libA_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libB_scp_meta_2, file="libB_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libC_scp_meta_2, file="libC_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libD_scp_meta_2, file="libD_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libE_scp_meta_2, file="libE_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libF_scp_meta_2, file="libF_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libG_scp_meta_2, file="libG_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libH_scp_meta_2, file="libH_scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)


##CHECK WHY ONE ROW MISSING IN META_2 (missing?) --> check this)

a <-libH_scp_meta$NAME
b <- libH_scp_meta_2$NAME

setdiff(a,b)


##libA "AAACCCAAGCGTTAGG"
##libB "AAACGAAAGGTCCCTG"
##libC "AAACCCAAGAAACTAC"
##libD  "AAACCCAAGGCGAAGG"
##libE  "AAACCCACAGACAAAT"
##libF  equal
##libG  "AAACGAACACCAGCCA"
##libH  "AAACCCACAGCTGAAG"

##Add lib identifier
libA_scp_meta_2$NAME <- paste(libA_scp_meta_2$NAME,"A",sep="-")
libB_scp_meta_2$NAME <- paste(libB_scp_meta_2$NAME,"B",sep="-")
libC_scp_meta_2$NAME <- paste(libC_scp_meta_2$NAME,"C",sep="-")
libD_scp_meta_2$NAME <- paste(libD_scp_meta_2$NAME,"D",sep="-")
libE_scp_meta_2$NAME <- paste(libE_scp_meta_2$NAME,"E",sep="-")
libF_scp_meta_2$NAME <- paste(libF_scp_meta_2$NAME,"F",sep="-")
libG_scp_meta_2$NAME <- paste(libG_scp_meta_2$NAME,"G",sep="-")
libH_scp_meta_2$NAME <- paste(libH_scp_meta_2$NAME,"H",sep="-")

##merge all data into one file

scp_meta <- rbind(libA_scp_meta_2, libB_scp_meta_2, libC_scp_meta_2, libD_scp_meta_2,
                  libE_scp_meta_2, libF_scp_meta_2, libG_scp_meta_2, libH_scp_meta_2)

scp_meta_2 <- scp_meta

scp_meta_2$Cell_Type <- gsub("-", "_", scp_meta$Cell_Type)
scp_meta_2$Cell_Type <- gsub("/", "_", scp_meta_2$Cell_Type)

##save data

#write.table(scp_meta, file="scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(scp_meta_2, file="scp_meta_2.txt", sep="\t", quote=FALSE, row.names=FALSE)


##check if original scp_meta_2 is equal = yes, so the original
#scp_out_tex.txt.gz matrix files as well as the sum_scp_out_tex.txt.gz files should work? Nothing works anymore!!

#identical(libA_scp_meta_2$NAME,libA_scp_meta_2_og$NAME)
#identical(libA_scp_meta_2$Count,libA_scp_meta_2_og$Count)

#identical(libB_scp_meta_2$NAME,libB_scp_meta_2_og$NAME)
#identical(libB_scp_meta_2$Count,libB_scp_meta_2_og$Count)

#identical(libC_scp_meta_2$NAME,libC_scp_meta_2_og$NAME)
#identical(libC_scp_meta_2$Count,libC_scp_meta_2_og$Count)

#identical(libD_scp_meta_2$NAME,libD_scp_meta_2_og$NAME)
#identical(libD_scp_meta_2$Count,libD_scp_meta_2_og$Count)

#identical(libE_scp_meta_2$NAME,libE_scp_meta_2_og$NAME)
#identical(libE_scp_meta_2$Count,libE_scp_meta_2_og$Count)

#identical(libF_scp_meta_2$NAME,libF_scp_meta_2_og$NAME)
#identical(libF_scp_meta_2$Count,libF_scp_meta_2_og$Count)

#identical(libG_scp_meta_2$NAME,libG_scp_meta_2_og$NAME)
#identical(libG_scp_meta_2$Count,libG_scp_meta_2_og$Count)

#identical(libH_scp_meta_2$NAME,libH_scp_meta_2_og$NAME)
#identical(libH_scp_meta_2$Count,libH_scp_meta_2_og$Count)


#identical(scp_meta$NAME, scp_meta_2$NAME)
#identical(scp_meta$Counts, scp_meta_2$Counts)
#identical(scp_meta$Cell_State, scp_meta_2$Cell_State)
#identical(scp_meta$Cohort, scp_meta_2$Cohort)
#identical(scp_meta$Patient, scp_meta_2$Patient)
