library(Seurat)
#library(ggplot2)
library(sctransform)
library(glmGamPoi)
#library(metap)
#library(multtest)
library(dplyr)
#library(data.table)
#library(stringr)
#library(dplyr)
#library(tidyverse)
#library(RCurl)
#library(cowplot)
#library(ggplot2)
#library(hrbrthemes)
#library(viridis)
library(dittoSeq)
#memory.limit(size=1600000)

#Check ALL DATASETS

Long_read_UK_72hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_UK_72hpi_adult/outs/filtered_feature_bc_matrix")
Long_read_UK_72hpi_adult <- CreateSeuratObject(counts = Long_read_UK_72hpi_adult)
demuxlet_Long_read_UK_72hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_UK_72hpi_adult/outs/demuxlet/Long_read_UK_72hpi_adult_nowarning.best", header=TRUE)


Long_read_UK_72hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_UK_72hpi_adult/raw_feature_bc_matrix")
Long_read_UK_72hpi_adult_v <- CreateSeuratObject(counts = Long_read_UK_72hpi_adult_v)


#Select infected cells
Long_read_UK_72hpi_adult_v@meta.data <- subset(Long_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_UK_72hpi_adult_v_barcode <- rownames(Long_read_UK_72hpi_adult_v@meta.data)
Long_read_UK_72hpi_adult_v_barcode <- as.list(Long_read_UK_72hpi_adult_v_barcode)

#Add infection level
Long_read_UK_72hpi_adult@meta.data$Infection <- ifelse(rownames(Long_read_UK_72hpi_adult@meta.data) %in% Long_read_UK_72hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection == "Uninfected",]))


Long_read_UK_72hpi_adult_low_meta.data <- subset(Long_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_UK_72hpi_adult_medium_meta.data <- subset(Long_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_UK_72hpi_adult_high_meta.data <- subset(Long_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_UK_72hpi_adult_very_high_meta.data <- subset(Long_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_UK_72hpi_adult_low_barcode <- rownames(Long_read_UK_72hpi_adult_low_meta.data)
Long_read_UK_72hpi_adult_low_barcode <- as.list(Long_read_UK_72hpi_adult_low_barcode)

Long_read_UK_72hpi_adult_medium_barcode <- rownames(Long_read_UK_72hpi_adult_medium_meta.data)
Long_read_UK_72hpi_adult_medium_barcode <- as.list(Long_read_UK_72hpi_adult_medium_barcode)

Long_read_UK_72hpi_adult_high_barcode <- rownames(Long_read_UK_72hpi_adult_high_meta.data)
Long_read_UK_72hpi_adult_high_barcode <- as.list(Long_read_UK_72hpi_adult_high_barcode)

Long_read_UK_72hpi_adult_very_high_barcode <- rownames(Long_read_UK_72hpi_adult_very_high_meta.data)
Long_read_UK_72hpi_adult_very_high_barcode <- as.list(Long_read_UK_72hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_adult@meta.data) %in% Long_read_UK_72hpi_adult_low_barcode, print("Low"), print("Unknown"))
Long_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_adult@meta.data) %in% Long_read_UK_72hpi_adult_medium_barcode & Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_UK_72hpi_adult@meta.data$Infection_tier)
Long_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_adult@meta.data) %in% Long_read_UK_72hpi_adult_high_barcode & Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Long_read_UK_72hpi_adult@meta.data$Infection_tier)
Long_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_adult@meta.data) %in% Long_read_UK_72hpi_adult_very_high_barcode & Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_UK_72hpi_adult@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_UK_72hpi_adult@meta.data[Long_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown",]))


Long_read_UK_72hpi_adult <-importDemux(
    Long_read_UK_72hpi_adult,
    demuxlet.best = demuxlet_Long_read_UK_72hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_UK_72hpi_adult)
demux.calls.summary(Long_read_UK_72hpi_adult)



Long_read_uninfected_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_uninfected_adult/outs/filtered_feature_bc_matrix")
Long_read_uninfected_adult <- CreateSeuratObject(counts = Long_read_uninfected_adult)
demuxlet_Long_read_uninfected_adult <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_uninfected_adult/outs/demuxlet/Long_read_uninfected_adult_nowarning.best", header=TRUE)



Long_read_uninfected_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_uninfected_adult/raw_feature_bc_matrix")
Long_read_uninfected_adult_v <- CreateSeuratObject(counts = Long_read_uninfected_adult_v)


#Select infected cells
Long_read_uninfected_adult_v@meta.data <- subset(Long_read_uninfected_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_uninfected_adult_v_barcode <- rownames(Long_read_uninfected_adult_v@meta.data)
Long_read_uninfected_adult_v_barcode <- as.list(Long_read_uninfected_adult_v_barcode)

#Add infection level
Long_read_uninfected_adult@meta.data$Infection <- ifelse(rownames(Long_read_uninfected_adult@meta.data) %in% Long_read_uninfected_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection == "Infected",]))
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection == "Uninfected",]))

Long_read_uninfected_adult_low_meta.data <- subset(Long_read_uninfected_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_uninfected_adult_medium_meta.data <- subset(Long_read_uninfected_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_uninfected_adult_high_meta.data <- subset(Long_read_uninfected_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_uninfected_adult_very_high_meta.data <- subset(Long_read_uninfected_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_uninfected_adult_low_barcode <- rownames(Long_read_uninfected_adult_low_meta.data)
Long_read_uninfected_adult_low_barcode <- as.list(Long_read_uninfected_adult_low_barcode)

Long_read_uninfected_adult_medium_barcode <- rownames(Long_read_uninfected_adult_medium_meta.data)
Long_read_uninfected_adult_medium_barcode <- as.list(Long_read_uninfected_adult_medium_barcode)

Long_read_uninfected_adult_high_barcode <- rownames(Long_read_uninfected_adult_high_meta.data)
Long_read_uninfected_adult_high_barcode <- as.list(Long_read_uninfected_adult_high_barcode)

Long_read_uninfected_adult_very_high_barcode <- rownames(Long_read_uninfected_adult_very_high_meta.data)
Long_read_uninfected_adult_very_high_barcode <- as.list(Long_read_uninfected_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_adult@meta.data) %in% Long_read_uninfected_adult_low_barcode, print("Low"), print("Unknown"))
Long_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_adult@meta.data) %in% Long_read_uninfected_adult_medium_barcode & Long_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_uninfected_adult@meta.data$Infection_tier)
Long_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_adult@meta.data) %in% Long_read_uninfected_adult_high_barcode & Long_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("High"),Long_read_uninfected_adult@meta.data$Infection_tier)
Long_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_adult@meta.data) %in% Long_read_uninfected_adult_very_high_barcode & Long_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_uninfected_adult@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_uninfected_adult@meta.data[Long_read_uninfected_adult@meta.data$Infection_tier == "Unknown",]))


Long_read_uninfected_adult <-importDemux(
    Long_read_uninfected_adult,
    demuxlet.best = demuxlet_Long_read_uninfected_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_uninfected_adult)
demux.calls.summary(Long_read_uninfected_adult)

Long_read_VIC01_48hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_VIC01_48hpi_adult/outs/filtered_feature_bc_matrix")
Long_read_VIC01_48hpi_adult <- CreateSeuratObject(counts = Long_read_VIC01_48hpi_adult)
demuxlet_Long_read_VIC01_48hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_VIC01_48hpi_adult/outs/demuxlet/Long_read_VIC01_48hpi_adult_nowarning.best", header=TRUE)



Long_read_VIC01_48hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_VIC01_48hpi_adult/raw_feature_bc_matrix")
Long_read_VIC01_48hpi_adult_v <- CreateSeuratObject(counts = Long_read_VIC01_48hpi_adult_v)


#Select infected cells
Long_read_VIC01_48hpi_adult_v@meta.data <- subset(Long_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_VIC01_48hpi_adult_v_barcode <- rownames(Long_read_VIC01_48hpi_adult_v@meta.data)
Long_read_VIC01_48hpi_adult_v_barcode <- as.list(Long_read_VIC01_48hpi_adult_v_barcode)

#Add infection level
Long_read_VIC01_48hpi_adult@meta.data$Infection <- ifelse(rownames(Long_read_VIC01_48hpi_adult@meta.data) %in% Long_read_VIC01_48hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection == "Uninfected",]))


Long_read_VIC01_48hpi_adult_low_meta.data <- subset(Long_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_VIC01_48hpi_adult_medium_meta.data <- subset(Long_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_VIC01_48hpi_adult_high_meta.data <- subset(Long_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_VIC01_48hpi_adult_very_high_meta.data <- subset(Long_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_VIC01_48hpi_adult_low_barcode <- rownames(Long_read_VIC01_48hpi_adult_low_meta.data)
Long_read_VIC01_48hpi_adult_low_barcode <- as.list(Long_read_VIC01_48hpi_adult_low_barcode)

Long_read_VIC01_48hpi_adult_medium_barcode <- rownames(Long_read_VIC01_48hpi_adult_medium_meta.data)
Long_read_VIC01_48hpi_adult_medium_barcode <- as.list(Long_read_VIC01_48hpi_adult_medium_barcode)

Long_read_VIC01_48hpi_adult_high_barcode <- rownames(Long_read_VIC01_48hpi_adult_high_meta.data)
Long_read_VIC01_48hpi_adult_high_barcode <- as.list(Long_read_VIC01_48hpi_adult_high_barcode)

Long_read_VIC01_48hpi_adult_very_high_barcode <- rownames(Long_read_VIC01_48hpi_adult_very_high_meta.data)
Long_read_VIC01_48hpi_adult_very_high_barcode <- as.list(Long_read_VIC01_48hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_adult@meta.data) %in% Long_read_VIC01_48hpi_adult_low_barcode, print("Low"), print("Unknown"))
Long_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_adult@meta.data) %in% Long_read_VIC01_48hpi_adult_medium_barcode & Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_VIC01_48hpi_adult@meta.data$Infection_tier)
Long_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_adult@meta.data) %in% Long_read_VIC01_48hpi_adult_high_barcode & Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Long_read_VIC01_48hpi_adult@meta.data$Infection_tier)
Long_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_adult@meta.data) %in% Long_read_VIC01_48hpi_adult_very_high_barcode & Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_VIC01_48hpi_adult@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_VIC01_48hpi_adult@meta.data[Long_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown",]))


Long_read_VIC01_48hpi_adult <-importDemux(
    Long_read_VIC01_48hpi_adult,
    demuxlet.best = demuxlet_Long_read_VIC01_48hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_VIC01_48hpi_adult)
demux.calls.summary(Long_read_VIC01_48hpi_adult)

Long_read_VIC01_72hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_VIC01_72hpi_adult/outs/filtered_feature_bc_matrix")
Long_read_VIC01_72hpi_adult <- CreateSeuratObject(counts = Long_read_VIC01_72hpi_adult)
demuxlet_Long_read_VIC01_72hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_VIC01_72hpi_adult/outs/demuxlet/Long_read_VIC01_72hpi_adult_nowarning.best", header=TRUE)


Long_read_VIC01_72hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_VIC01_72hpi_adult/raw_feature_bc_matrix")
Long_read_VIC01_72hpi_adult_v <- CreateSeuratObject(counts = Long_read_VIC01_72hpi_adult_v)


#Select infected cells
Long_read_VIC01_72hpi_adult_v@meta.data <- subset(Long_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_VIC01_72hpi_adult_v_barcode <- rownames(Long_read_VIC01_72hpi_adult_v@meta.data)
Long_read_VIC01_72hpi_adult_v_barcode <- as.list(Long_read_VIC01_72hpi_adult_v_barcode)

#Add infection level
Long_read_VIC01_72hpi_adult@meta.data$Infection <- ifelse(rownames(Long_read_VIC01_72hpi_adult@meta.data) %in% Long_read_VIC01_72hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection == "Uninfected",]))


Long_read_VIC01_72hpi_adult_low_meta.data <- subset(Long_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_VIC01_72hpi_adult_medium_meta.data <- subset(Long_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_VIC01_72hpi_adult_high_meta.data <- subset(Long_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_VIC01_72hpi_adult_very_high_meta.data <- subset(Long_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_VIC01_72hpi_adult_low_barcode <- rownames(Long_read_VIC01_72hpi_adult_low_meta.data)
Long_read_VIC01_72hpi_adult_low_barcode <- as.list(Long_read_VIC01_72hpi_adult_low_barcode)

Long_read_VIC01_72hpi_adult_medium_barcode <- rownames(Long_read_VIC01_72hpi_adult_medium_meta.data)
Long_read_VIC01_72hpi_adult_medium_barcode <- as.list(Long_read_VIC01_72hpi_adult_medium_barcode)

Long_read_VIC01_72hpi_adult_high_barcode <- rownames(Long_read_VIC01_72hpi_adult_high_meta.data)
Long_read_VIC01_72hpi_adult_high_barcode <- as.list(Long_read_VIC01_72hpi_adult_high_barcode)

Long_read_VIC01_72hpi_adult_very_high_barcode <- rownames(Long_read_VIC01_72hpi_adult_very_high_meta.data)
Long_read_VIC01_72hpi_adult_very_high_barcode <- as.list(Long_read_VIC01_72hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_adult@meta.data) %in% Long_read_VIC01_72hpi_adult_low_barcode, print("Low"), print("Unknown"))
Long_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_adult@meta.data) %in% Long_read_VIC01_72hpi_adult_medium_barcode & Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_VIC01_72hpi_adult@meta.data$Infection_tier)
Long_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_adult@meta.data) %in% Long_read_VIC01_72hpi_adult_high_barcode & Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Long_read_VIC01_72hpi_adult@meta.data$Infection_tier)
Long_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_adult@meta.data) %in% Long_read_VIC01_72hpi_adult_very_high_barcode & Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_VIC01_72hpi_adult@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_VIC01_72hpi_adult@meta.data[Long_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown",]))


Long_read_VIC01_72hpi_adult <-importDemux(
    Long_read_VIC01_72hpi_adult,
    demuxlet.best = demuxlet_Long_read_VIC01_72hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Long_read_VIC01_72hpi_adult)
demux.calls.summary(Long_read_VIC01_72hpi_adult)


Short_read_UK_72hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_UK_72hpi_adult/outs/filtered_feature_bc_matrix")
Short_read_UK_72hpi_adult <- CreateSeuratObject(counts = Short_read_UK_72hpi_adult)
demuxlet_Short_read_UK_72hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_UK_72hpi_adult/outs/demuxlet/Short_read_UK_72hpi_adult_nowarning.best", header=TRUE)



Short_read_UK_72hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_UK_72hpi_adult/raw_feature_bc_matrix")
Short_read_UK_72hpi_adult_v <- CreateSeuratObject(counts = Short_read_UK_72hpi_adult_v)


#Select infected cells
Short_read_UK_72hpi_adult_v@meta.data <- subset(Short_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_UK_72hpi_adult_v_barcode <- rownames(Short_read_UK_72hpi_adult_v@meta.data)
Short_read_UK_72hpi_adult_v_barcode <- as.list(Short_read_UK_72hpi_adult_v_barcode)

#Add infection level
Short_read_UK_72hpi_adult@meta.data$Infection <- ifelse(rownames(Short_read_UK_72hpi_adult@meta.data) %in% Short_read_UK_72hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection == "Uninfected",]))

Short_read_UK_72hpi_adult_low_meta.data <- subset(Short_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_UK_72hpi_adult_medium_meta.data <- subset(Short_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_UK_72hpi_adult_high_meta.data <- subset(Short_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_UK_72hpi_adult_very_high_meta.data <- subset(Short_read_UK_72hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_UK_72hpi_adult_low_barcode <- rownames(Short_read_UK_72hpi_adult_low_meta.data)
Short_read_UK_72hpi_adult_low_barcode <- as.list(Short_read_UK_72hpi_adult_low_barcode)

Short_read_UK_72hpi_adult_medium_barcode <- rownames(Short_read_UK_72hpi_adult_medium_meta.data)
Short_read_UK_72hpi_adult_medium_barcode <- as.list(Short_read_UK_72hpi_adult_medium_barcode)

Short_read_UK_72hpi_adult_high_barcode <- rownames(Short_read_UK_72hpi_adult_high_meta.data)
Short_read_UK_72hpi_adult_high_barcode <- as.list(Short_read_UK_72hpi_adult_high_barcode)

Short_read_UK_72hpi_adult_very_high_barcode <- rownames(Short_read_UK_72hpi_adult_very_high_meta.data)
Short_read_UK_72hpi_adult_very_high_barcode <- as.list(Short_read_UK_72hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_adult@meta.data) %in% Short_read_UK_72hpi_adult_low_barcode, print("Low"), print("Unknown"))
Short_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_adult@meta.data) %in% Short_read_UK_72hpi_adult_medium_barcode & Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_UK_72hpi_adult@meta.data$Infection_tier)
Short_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_adult@meta.data) %in% Short_read_UK_72hpi_adult_high_barcode & Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Short_read_UK_72hpi_adult@meta.data$Infection_tier)
Short_read_UK_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_adult@meta.data) %in% Short_read_UK_72hpi_adult_very_high_barcode & Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_UK_72hpi_adult@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_UK_72hpi_adult@meta.data[Short_read_UK_72hpi_adult@meta.data$Infection_tier == "Unknown",]))



Short_read_UK_72hpi_adult <-importDemux(
    Short_read_UK_72hpi_adult,
    demuxlet.best = demuxlet_Short_read_UK_72hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Short_read_UK_72hpi_adult)
demux.calls.summary(Short_read_UK_72hpi_adult)


Short_read_uninfected_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_uninfected_adult/outs/filtered_feature_bc_matrix")
Short_read_uninfected_adult <- CreateSeuratObject(counts = Short_read_uninfected_adult)
demuxlet_Short_read_uninfected_adult <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_uninfected_adult/outs/demuxlet/Short_read_uninfected_adult_nowarning.best", header=TRUE)



Short_read_uninfected_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_uninfected_adult/raw_feature_bc_matrix")
Short_read_uninfected_adult_v <- CreateSeuratObject(counts = Short_read_uninfected_adult_v)


#Select infected cells
Short_read_uninfected_adult_v@meta.data <- subset(Short_read_uninfected_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_uninfected_adult_v_barcode <- rownames(Short_read_uninfected_adult_v@meta.data)
Short_read_uninfected_adult_v_barcode <- as.list(Short_read_uninfected_adult_v_barcode)

#Add infection level
Short_read_uninfected_adult@meta.data$Infection <- ifelse(rownames(Short_read_uninfected_adult@meta.data) %in% Short_read_uninfected_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection == "Infected",]))
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection == "Uninfected",]))


Short_read_uninfected_adult_low_meta.data <- subset(Short_read_uninfected_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_uninfected_adult_medium_meta.data <- subset(Short_read_uninfected_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_uninfected_adult_high_meta.data <- subset(Short_read_uninfected_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_uninfected_adult_very_high_meta.data <- subset(Short_read_uninfected_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_uninfected_adult_low_barcode <- rownames(Short_read_uninfected_adult_low_meta.data)
Short_read_uninfected_adult_low_barcode <- as.list(Short_read_uninfected_adult_low_barcode)

Short_read_uninfected_adult_medium_barcode <- rownames(Short_read_uninfected_adult_medium_meta.data)
Short_read_uninfected_adult_medium_barcode <- as.list(Short_read_uninfected_adult_medium_barcode)

Short_read_uninfected_adult_high_barcode <- rownames(Short_read_uninfected_adult_high_meta.data)
Short_read_uninfected_adult_high_barcode <- as.list(Short_read_uninfected_adult_high_barcode)

Short_read_uninfected_adult_very_high_barcode <- rownames(Short_read_uninfected_adult_very_high_meta.data)
Short_read_uninfected_adult_very_high_barcode <- as.list(Short_read_uninfected_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_adult@meta.data) %in% Short_read_uninfected_adult_low_barcode, print("Low"), print("Unknown"))
Short_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_adult@meta.data) %in% Short_read_uninfected_adult_medium_barcode & Short_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_uninfected_adult@meta.data$Infection_tier)
Short_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_adult@meta.data) %in% Short_read_uninfected_adult_high_barcode & Short_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("High"),Short_read_uninfected_adult@meta.data$Infection_tier)
Short_read_uninfected_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_adult@meta.data) %in% Short_read_uninfected_adult_very_high_barcode & Short_read_uninfected_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_uninfected_adult@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_uninfected_adult@meta.data[Short_read_uninfected_adult@meta.data$Infection_tier == "Unknown",]))



Short_read_uninfected_adult <-importDemux(
    Short_read_uninfected_adult,
    demuxlet.best = demuxlet_Short_read_uninfected_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Short_read_uninfected_adult)
demux.calls.summary(Short_read_uninfected_adult)



Short_read_VIC01_48hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_VIC01_48hpi_adult/outs/filtered_feature_bc_matrix")
Short_read_VIC01_48hpi_adult <- CreateSeuratObject(counts = Short_read_VIC01_48hpi_adult)
demuxlet_Short_read_VIC01_48hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_VIC01_48hpi_adult/outs/demuxlet/Short_read_VIC01_48hpi_adult_nowarning.best", header=TRUE)



Short_read_VIC01_48hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_VIC01_48hpi_adult/raw_feature_bc_matrix")
Short_read_VIC01_48hpi_adult_v <- CreateSeuratObject(counts = Short_read_VIC01_48hpi_adult_v)


#Select infected cells
Short_read_VIC01_48hpi_adult_v@meta.data <- subset(Short_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_VIC01_48hpi_adult_v_barcode <- rownames(Short_read_VIC01_48hpi_adult_v@meta.data)
Short_read_VIC01_48hpi_adult_v_barcode <- as.list(Short_read_VIC01_48hpi_adult_v_barcode)

#Add infection level
Short_read_VIC01_48hpi_adult@meta.data$Infection <- ifelse(rownames(Short_read_VIC01_48hpi_adult@meta.data) %in% Short_read_VIC01_48hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection == "Uninfected",]))

Short_read_VIC01_48hpi_adult_low_meta.data <- subset(Short_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_VIC01_48hpi_adult_medium_meta.data <- subset(Short_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_VIC01_48hpi_adult_high_meta.data <- subset(Short_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_VIC01_48hpi_adult_very_high_meta.data <- subset(Short_read_VIC01_48hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_VIC01_48hpi_adult_low_barcode <- rownames(Short_read_VIC01_48hpi_adult_low_meta.data)
Short_read_VIC01_48hpi_adult_low_barcode <- as.list(Short_read_VIC01_48hpi_adult_low_barcode)

Short_read_VIC01_48hpi_adult_medium_barcode <- rownames(Short_read_VIC01_48hpi_adult_medium_meta.data)
Short_read_VIC01_48hpi_adult_medium_barcode <- as.list(Short_read_VIC01_48hpi_adult_medium_barcode)

Short_read_VIC01_48hpi_adult_high_barcode <- rownames(Short_read_VIC01_48hpi_adult_high_meta.data)
Short_read_VIC01_48hpi_adult_high_barcode <- as.list(Short_read_VIC01_48hpi_adult_high_barcode)

Short_read_VIC01_48hpi_adult_very_high_barcode <- rownames(Short_read_VIC01_48hpi_adult_very_high_meta.data)
Short_read_VIC01_48hpi_adult_very_high_barcode <- as.list(Short_read_VIC01_48hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_adult@meta.data) %in% Short_read_VIC01_48hpi_adult_low_barcode, print("Low"), print("Unknown"))
Short_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_adult@meta.data) %in% Short_read_VIC01_48hpi_adult_medium_barcode & Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_VIC01_48hpi_adult@meta.data$Infection_tier)
Short_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_adult@meta.data) %in% Short_read_VIC01_48hpi_adult_high_barcode & Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Short_read_VIC01_48hpi_adult@meta.data$Infection_tier)
Short_read_VIC01_48hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_adult@meta.data) %in% Short_read_VIC01_48hpi_adult_very_high_barcode & Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_VIC01_48hpi_adult@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_VIC01_48hpi_adult@meta.data[Short_read_VIC01_48hpi_adult@meta.data$Infection_tier == "Unknown",]))



Short_read_VIC01_48hpi_adult <-importDemux(
    Short_read_VIC01_48hpi_adult,
    demuxlet.best = demuxlet_Short_read_VIC01_48hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Short_read_VIC01_48hpi_adult)
demux.calls.summary(Short_read_VIC01_48hpi_adult)



Short_read_VIC01_72hpi_adult <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_VIC01_72hpi_adult/outs/filtered_feature_bc_matrix")
Short_read_VIC01_72hpi_adult <- CreateSeuratObject(counts = Short_read_VIC01_72hpi_adult)
demuxlet_Short_read_VIC01_72hpi_adult <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_VIC01_72hpi_adult/outs/demuxlet/Short_read_VIC01_72hpi_adult_nowarning.best", header=TRUE)


Short_read_VIC01_72hpi_adult_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_VIC01_72hpi_adult/raw_feature_bc_matrix")
Short_read_VIC01_72hpi_adult_v <- CreateSeuratObject(counts = Short_read_VIC01_72hpi_adult_v)


#Select infected cells
Short_read_VIC01_72hpi_adult_v@meta.data <- subset(Short_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_VIC01_72hpi_adult_v_barcode <- rownames(Short_read_VIC01_72hpi_adult_v@meta.data)
Short_read_VIC01_72hpi_adult_v_barcode <- as.list(Short_read_VIC01_72hpi_adult_v_barcode)

#Add infection level
Short_read_VIC01_72hpi_adult@meta.data$Infection <- ifelse(rownames(Short_read_VIC01_72hpi_adult@meta.data) %in% Short_read_VIC01_72hpi_adult_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection == "Infected",]))
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection == "Uninfected",]))

Short_read_VIC01_72hpi_adult_low_meta.data <- subset(Short_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_VIC01_72hpi_adult_medium_meta.data <- subset(Short_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_VIC01_72hpi_adult_high_meta.data <- subset(Short_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_VIC01_72hpi_adult_very_high_meta.data <- subset(Short_read_VIC01_72hpi_adult_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_VIC01_72hpi_adult_low_barcode <- rownames(Short_read_VIC01_72hpi_adult_low_meta.data)
Short_read_VIC01_72hpi_adult_low_barcode <- as.list(Short_read_VIC01_72hpi_adult_low_barcode)

Short_read_VIC01_72hpi_adult_medium_barcode <- rownames(Short_read_VIC01_72hpi_adult_medium_meta.data)
Short_read_VIC01_72hpi_adult_medium_barcode <- as.list(Short_read_VIC01_72hpi_adult_medium_barcode)

Short_read_VIC01_72hpi_adult_high_barcode <- rownames(Short_read_VIC01_72hpi_adult_high_meta.data)
Short_read_VIC01_72hpi_adult_high_barcode <- as.list(Short_read_VIC01_72hpi_adult_high_barcode)

Short_read_VIC01_72hpi_adult_very_high_barcode <- rownames(Short_read_VIC01_72hpi_adult_very_high_meta.data)
Short_read_VIC01_72hpi_adult_very_high_barcode <- as.list(Short_read_VIC01_72hpi_adult_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_adult@meta.data) %in% Short_read_VIC01_72hpi_adult_low_barcode, print("Low"), print("Unknown"))
Short_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_adult@meta.data) %in% Short_read_VIC01_72hpi_adult_medium_barcode & Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_VIC01_72hpi_adult@meta.data$Infection_tier)
Short_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_adult@meta.data) %in% Short_read_VIC01_72hpi_adult_high_barcode & Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("High"),Short_read_VIC01_72hpi_adult@meta.data$Infection_tier)
Short_read_VIC01_72hpi_adult@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_adult@meta.data) %in% Short_read_VIC01_72hpi_adult_very_high_barcode & Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_VIC01_72hpi_adult@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_VIC01_72hpi_adult@meta.data[Short_read_VIC01_72hpi_adult@meta.data$Infection_tier == "Unknown",]))



Short_read_VIC01_72hpi_adult <-importDemux(
    Short_read_VIC01_72hpi_adult,
    demuxlet.best = demuxlet_Short_read_VIC01_72hpi_adult,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Short_read_VIC01_72hpi_adult)
demux.calls.summary(Short_read_VIC01_72hpi_adult)



Long_read_UK_72hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_UK_72hpi_child/outs/filtered_feature_bc_matrix")
Long_read_UK_72hpi_child <- CreateSeuratObject(counts = Long_read_UK_72hpi_child)
demuxlet_Long_read_UK_72hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_UK_72hpi_child/outs/demuxlet/Long_read_UK_72hpi_child_nowarning.best", header=TRUE)


Long_read_UK_72hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_UK_72hpi_child/raw_feature_bc_matrix")
Long_read_UK_72hpi_child_v <- CreateSeuratObject(counts = Long_read_UK_72hpi_child_v)


#Select infected cells
Long_read_UK_72hpi_child_v@meta.data <- subset(Long_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_UK_72hpi_child_v_barcode <- rownames(Long_read_UK_72hpi_child_v@meta.data)
Long_read_UK_72hpi_child_v_barcode <- as.list(Long_read_UK_72hpi_child_v_barcode)

#Add infection level
Long_read_UK_72hpi_child@meta.data$Infection <- ifelse(rownames(Long_read_UK_72hpi_child@meta.data) %in% Long_read_UK_72hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection == "Uninfected",]))

Long_read_UK_72hpi_child_low_meta.data <- subset(Long_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_UK_72hpi_child_medium_meta.data <- subset(Long_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_UK_72hpi_child_high_meta.data <- subset(Long_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_UK_72hpi_child_very_high_meta.data <- subset(Long_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_UK_72hpi_child_low_barcode <- rownames(Long_read_UK_72hpi_child_low_meta.data)
Long_read_UK_72hpi_child_low_barcode <- as.list(Long_read_UK_72hpi_child_low_barcode)

Long_read_UK_72hpi_child_medium_barcode <- rownames(Long_read_UK_72hpi_child_medium_meta.data)
Long_read_UK_72hpi_child_medium_barcode <- as.list(Long_read_UK_72hpi_child_medium_barcode)

Long_read_UK_72hpi_child_high_barcode <- rownames(Long_read_UK_72hpi_child_high_meta.data)
Long_read_UK_72hpi_child_high_barcode <- as.list(Long_read_UK_72hpi_child_high_barcode)

Long_read_UK_72hpi_child_very_high_barcode <- rownames(Long_read_UK_72hpi_child_very_high_meta.data)
Long_read_UK_72hpi_child_very_high_barcode <- as.list(Long_read_UK_72hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_child@meta.data) %in% Long_read_UK_72hpi_child_low_barcode, print("Low"), print("Unknown"))
Long_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_child@meta.data) %in% Long_read_UK_72hpi_child_medium_barcode & Long_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_UK_72hpi_child@meta.data$Infection_tier)
Long_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_child@meta.data) %in% Long_read_UK_72hpi_child_high_barcode & Long_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Long_read_UK_72hpi_child@meta.data$Infection_tier)
Long_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_UK_72hpi_child@meta.data) %in% Long_read_UK_72hpi_child_very_high_barcode & Long_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_UK_72hpi_child@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_UK_72hpi_child@meta.data[Long_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown",]))


Long_read_UK_72hpi_child <-importDemux(
    Long_read_UK_72hpi_child,
    demuxlet.best = demuxlet_Long_read_UK_72hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_UK_72hpi_child)
demux.calls.summary(Long_read_UK_72hpi_child)



Long_read_uninfected_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_uninfected_child/outs/filtered_feature_bc_matrix")
Long_read_uninfected_child <- CreateSeuratObject(counts = Long_read_uninfected_child)
demuxlet_Long_read_uninfected_child <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_uninfected_child/outs/demuxlet/Long_read_uninfected_child_nowarning.best", header=TRUE)

Long_read_uninfected_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_uninfected_child/raw_feature_bc_matrix")
Long_read_uninfected_child_v <- CreateSeuratObject(counts = Long_read_uninfected_child_v)


#Select infected cells
Long_read_uninfected_child_v@meta.data <- subset(Long_read_uninfected_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_uninfected_child_v_barcode <- rownames(Long_read_uninfected_child_v@meta.data)
Long_read_uninfected_child_v_barcode <- as.list(Long_read_uninfected_child_v_barcode)

#Add infection level
Long_read_uninfected_child@meta.data$Infection <- ifelse(rownames(Long_read_uninfected_child@meta.data) %in% Long_read_uninfected_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection == "Infected",]))
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection == "Uninfected",]))

Long_read_uninfected_child_low_meta.data <- subset(Long_read_uninfected_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_uninfected_child_medium_meta.data <- subset(Long_read_uninfected_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_uninfected_child_high_meta.data <- subset(Long_read_uninfected_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_uninfected_child_very_high_meta.data <- subset(Long_read_uninfected_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_uninfected_child_low_barcode <- rownames(Long_read_uninfected_child_low_meta.data)
Long_read_uninfected_child_low_barcode <- as.list(Long_read_uninfected_child_low_barcode)

Long_read_uninfected_child_medium_barcode <- rownames(Long_read_uninfected_child_medium_meta.data)
Long_read_uninfected_child_medium_barcode <- as.list(Long_read_uninfected_child_medium_barcode)

Long_read_uninfected_child_high_barcode <- rownames(Long_read_uninfected_child_high_meta.data)
Long_read_uninfected_child_high_barcode <- as.list(Long_read_uninfected_child_high_barcode)

Long_read_uninfected_child_very_high_barcode <- rownames(Long_read_uninfected_child_very_high_meta.data)
Long_read_uninfected_child_very_high_barcode <- as.list(Long_read_uninfected_child_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_child@meta.data) %in% Long_read_uninfected_child_low_barcode, print("Low"), print("Unknown"))
Long_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_child@meta.data) %in% Long_read_uninfected_child_medium_barcode & Long_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_uninfected_child@meta.data$Infection_tier)
Long_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_child@meta.data) %in% Long_read_uninfected_child_high_barcode & Long_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("High"),Long_read_uninfected_child@meta.data$Infection_tier)
Long_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_uninfected_child@meta.data) %in% Long_read_uninfected_child_very_high_barcode & Long_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_uninfected_child@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_uninfected_child@meta.data[Long_read_uninfected_child@meta.data$Infection_tier == "Unknown",]))


Long_read_uninfected_child <-importDemux(
    Long_read_uninfected_child,
    demuxlet.best = demuxlet_Long_read_uninfected_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_uninfected_child)
demux.calls.summary(Long_read_uninfected_child)


Long_read_VIC01_48hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_VIC01_48hpi_child/outs/filtered_feature_bc_matrix")
Long_read_VIC01_48hpi_child <- CreateSeuratObject(counts = Long_read_VIC01_48hpi_child)
demuxlet_Long_read_VIC01_48hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_VIC01_48hpi_child/outs/demuxlet/Long_read_VIC01_48hpi_child_nowarning.best", header=TRUE)

Long_read_VIC01_48hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_VIC01_48hpi_child/raw_feature_bc_matrix")
Long_read_VIC01_48hpi_child_v <- CreateSeuratObject(counts = Long_read_VIC01_48hpi_child_v)


#Select infected cells
Long_read_VIC01_48hpi_child_v@meta.data <- subset(Long_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_VIC01_48hpi_child_v_barcode <- rownames(Long_read_VIC01_48hpi_child_v@meta.data)
Long_read_VIC01_48hpi_child_v_barcode <- as.list(Long_read_VIC01_48hpi_child_v_barcode)

#Add infection level
Long_read_VIC01_48hpi_child@meta.data$Infection <- ifelse(rownames(Long_read_VIC01_48hpi_child@meta.data) %in% Long_read_VIC01_48hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection == "Uninfected",]))

Long_read_VIC01_48hpi_child_low_meta.data <- subset(Long_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_VIC01_48hpi_child_medium_meta.data <- subset(Long_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_VIC01_48hpi_child_high_meta.data <- subset(Long_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_VIC01_48hpi_child_very_high_meta.data <- subset(Long_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_VIC01_48hpi_child_low_barcode <- rownames(Long_read_VIC01_48hpi_child_low_meta.data)
Long_read_VIC01_48hpi_child_low_barcode <- as.list(Long_read_VIC01_48hpi_child_low_barcode)

Long_read_VIC01_48hpi_child_medium_barcode <- rownames(Long_read_VIC01_48hpi_child_medium_meta.data)
Long_read_VIC01_48hpi_child_medium_barcode <- as.list(Long_read_VIC01_48hpi_child_medium_barcode)

Long_read_VIC01_48hpi_child_high_barcode <- rownames(Long_read_VIC01_48hpi_child_high_meta.data)
Long_read_VIC01_48hpi_child_high_barcode <- as.list(Long_read_VIC01_48hpi_child_high_barcode)

Long_read_VIC01_48hpi_child_very_high_barcode <- rownames(Long_read_VIC01_48hpi_child_very_high_meta.data)
Long_read_VIC01_48hpi_child_very_high_barcode <- as.list(Long_read_VIC01_48hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_child@meta.data) %in% Long_read_VIC01_48hpi_child_low_barcode, print("Low"), print("Unknown"))
Long_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_child@meta.data) %in% Long_read_VIC01_48hpi_child_medium_barcode & Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_VIC01_48hpi_child@meta.data$Infection_tier)
Long_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_child@meta.data) %in% Long_read_VIC01_48hpi_child_high_barcode & Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Long_read_VIC01_48hpi_child@meta.data$Infection_tier)
Long_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_48hpi_child@meta.data) %in% Long_read_VIC01_48hpi_child_very_high_barcode & Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_VIC01_48hpi_child@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_VIC01_48hpi_child@meta.data[Long_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown",]))



Long_read_VIC01_48hpi_child <-importDemux(
    Long_read_VIC01_48hpi_child,
    demuxlet.best = demuxlet_Long_read_VIC01_48hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Long_read_VIC01_48hpi_child)
demux.calls.summary(Long_read_VIC01_48hpi_child)



Long_read_VIC01_72hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Long_read_VIC01_72hpi_child/outs/filtered_feature_bc_matrix")
Long_read_VIC01_72hpi_child <- CreateSeuratObject(counts = Long_read_VIC01_72hpi_child)
demuxlet_Long_read_VIC01_72hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Long_read_VIC01_72hpi_child/outs/demuxlet/Long_read_VIC01_72hpi_child_nowarning.best", header=TRUE)


Long_read_VIC01_72hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Long_read_VIC01_72hpi_child/raw_feature_bc_matrix")
Long_read_VIC01_72hpi_child_v <- CreateSeuratObject(counts = Long_read_VIC01_72hpi_child_v)


#Select infected cells
Long_read_VIC01_72hpi_child_v@meta.data <- subset(Long_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Long_read_VIC01_72hpi_child_v_barcode <- rownames(Long_read_VIC01_72hpi_child_v@meta.data)
Long_read_VIC01_72hpi_child_v_barcode <- as.list(Long_read_VIC01_72hpi_child_v_barcode)

#Add infection level
Long_read_VIC01_72hpi_child@meta.data$Infection <- ifelse(rownames(Long_read_VIC01_72hpi_child@meta.data) %in% Long_read_VIC01_72hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection == "Uninfected",]))

Long_read_VIC01_72hpi_child_low_meta.data <- subset(Long_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Long_read_VIC01_72hpi_child_medium_meta.data <- subset(Long_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Long_read_VIC01_72hpi_child_high_meta.data <- subset(Long_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Long_read_VIC01_72hpi_child_very_high_meta.data <- subset(Long_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Long_read_VIC01_72hpi_child_low_barcode <- rownames(Long_read_VIC01_72hpi_child_low_meta.data)
Long_read_VIC01_72hpi_child_low_barcode <- as.list(Long_read_VIC01_72hpi_child_low_barcode)

Long_read_VIC01_72hpi_child_medium_barcode <- rownames(Long_read_VIC01_72hpi_child_medium_meta.data)
Long_read_VIC01_72hpi_child_medium_barcode <- as.list(Long_read_VIC01_72hpi_child_medium_barcode)

Long_read_VIC01_72hpi_child_high_barcode <- rownames(Long_read_VIC01_72hpi_child_high_meta.data)
Long_read_VIC01_72hpi_child_high_barcode <- as.list(Long_read_VIC01_72hpi_child_high_barcode)

Long_read_VIC01_72hpi_child_very_high_barcode <- rownames(Long_read_VIC01_72hpi_child_very_high_meta.data)
Long_read_VIC01_72hpi_child_very_high_barcode <- as.list(Long_read_VIC01_72hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Long_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_child@meta.data) %in% Long_read_VIC01_72hpi_child_low_barcode, print("Low"), print("Unknown"))
Long_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_child@meta.data) %in% Long_read_VIC01_72hpi_child_medium_barcode & Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Long_read_VIC01_72hpi_child@meta.data$Infection_tier)
Long_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_child@meta.data) %in% Long_read_VIC01_72hpi_child_high_barcode & Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Long_read_VIC01_72hpi_child@meta.data$Infection_tier)
Long_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Long_read_VIC01_72hpi_child@meta.data) %in% Long_read_VIC01_72hpi_child_very_high_barcode & Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Long_read_VIC01_72hpi_child@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Long_read_VIC01_72hpi_child@meta.data[Long_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown",]))



Long_read_VIC01_72hpi_child <-importDemux(
    Long_read_VIC01_72hpi_child,
    demuxlet.best = demuxlet_Long_read_VIC01_72hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Long_read_VIC01_72hpi_child)
demux.calls.summary(Long_read_VIC01_72hpi_child)


Short_read_UK_72hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_UK_72hpi_child/outs/filtered_feature_bc_matrix")
Short_read_UK_72hpi_child <- CreateSeuratObject(counts = Short_read_UK_72hpi_child)
demuxlet_Short_read_UK_72hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_UK_72hpi_child/outs/demuxlet/Short_read_UK_72hpi_child_nowarning.best", header=TRUE)


Short_read_UK_72hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_UK_72hpi_child/raw_feature_bc_matrix")
Short_read_UK_72hpi_child_v <- CreateSeuratObject(counts = Short_read_UK_72hpi_child_v)


#Select infected cells
Short_read_UK_72hpi_child_v@meta.data <- subset(Short_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_UK_72hpi_child_v_barcode <- rownames(Short_read_UK_72hpi_child_v@meta.data)
Short_read_UK_72hpi_child_v_barcode <- as.list(Short_read_UK_72hpi_child_v_barcode)

#Add infection level
Short_read_UK_72hpi_child@meta.data$Infection <- ifelse(rownames(Short_read_UK_72hpi_child@meta.data) %in% Short_read_UK_72hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection == "Uninfected",]))

Short_read_UK_72hpi_child_low_meta.data <- subset(Short_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_UK_72hpi_child_medium_meta.data <- subset(Short_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_UK_72hpi_child_high_meta.data <- subset(Short_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_UK_72hpi_child_very_high_meta.data <- subset(Short_read_UK_72hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_UK_72hpi_child_low_barcode <- rownames(Short_read_UK_72hpi_child_low_meta.data)
Short_read_UK_72hpi_child_low_barcode <- as.list(Short_read_UK_72hpi_child_low_barcode)

Short_read_UK_72hpi_child_medium_barcode <- rownames(Short_read_UK_72hpi_child_medium_meta.data)
Short_read_UK_72hpi_child_medium_barcode <- as.list(Short_read_UK_72hpi_child_medium_barcode)

Short_read_UK_72hpi_child_high_barcode <- rownames(Short_read_UK_72hpi_child_high_meta.data)
Short_read_UK_72hpi_child_high_barcode <- as.list(Short_read_UK_72hpi_child_high_barcode)

Short_read_UK_72hpi_child_very_high_barcode <- rownames(Short_read_UK_72hpi_child_very_high_meta.data)
Short_read_UK_72hpi_child_very_high_barcode <- as.list(Short_read_UK_72hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_child@meta.data) %in% Short_read_UK_72hpi_child_low_barcode, print("Low"), print("Unknown"))
Short_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_child@meta.data) %in% Short_read_UK_72hpi_child_medium_barcode & Short_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_UK_72hpi_child@meta.data$Infection_tier)
Short_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_child@meta.data) %in% Short_read_UK_72hpi_child_high_barcode & Short_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Short_read_UK_72hpi_child@meta.data$Infection_tier)
Short_read_UK_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_UK_72hpi_child@meta.data) %in% Short_read_UK_72hpi_child_very_high_barcode & Short_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_UK_72hpi_child@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_UK_72hpi_child@meta.data[Short_read_UK_72hpi_child@meta.data$Infection_tier == "Unknown",]))


Short_read_UK_72hpi_child <-importDemux(
    Short_read_UK_72hpi_child,
    demuxlet.best = demuxlet_Short_read_UK_72hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Short_read_UK_72hpi_child)
demux.calls.summary(Short_read_UK_72hpi_child)



Short_read_uninfected_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_uninfected_child/outs/filtered_feature_bc_matrix")
Short_read_uninfected_child <- CreateSeuratObject(counts = Short_read_uninfected_child)
demuxlet_Short_read_uninfected_child <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_uninfected_child/outs/demuxlet/Short_read_uninfected_child_nowarning.best", header=TRUE)

Short_read_uninfected_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_uninfected_child/raw_feature_bc_matrix")
Short_read_uninfected_child_v <- CreateSeuratObject(counts = Short_read_uninfected_child_v)


#Select infected cells
Short_read_uninfected_child_v@meta.data <- subset(Short_read_uninfected_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_uninfected_child_v_barcode <- rownames(Short_read_uninfected_child_v@meta.data)
Short_read_uninfected_child_v_barcode <- as.list(Short_read_uninfected_child_v_barcode)

#Add infection level
Short_read_uninfected_child@meta.data$Infection <- ifelse(rownames(Short_read_uninfected_child@meta.data) %in% Short_read_uninfected_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection == "Infected",]))
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection == "Uninfected",]))


Short_read_uninfected_child_low_meta.data <- subset(Short_read_uninfected_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_uninfected_child_medium_meta.data <- subset(Short_read_uninfected_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_uninfected_child_high_meta.data <- subset(Short_read_uninfected_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_uninfected_child_very_high_meta.data <- subset(Short_read_uninfected_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_uninfected_child_low_barcode <- rownames(Short_read_uninfected_child_low_meta.data)
Short_read_uninfected_child_low_barcode <- as.list(Short_read_uninfected_child_low_barcode)

Short_read_uninfected_child_medium_barcode <- rownames(Short_read_uninfected_child_medium_meta.data)
Short_read_uninfected_child_medium_barcode <- as.list(Short_read_uninfected_child_medium_barcode)

Short_read_uninfected_child_high_barcode <- rownames(Short_read_uninfected_child_high_meta.data)
Short_read_uninfected_child_high_barcode <- as.list(Short_read_uninfected_child_high_barcode)

Short_read_uninfected_child_very_high_barcode <- rownames(Short_read_uninfected_child_very_high_meta.data)
Short_read_uninfected_child_very_high_barcode <- as.list(Short_read_uninfected_child_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_child@meta.data) %in% Short_read_uninfected_child_low_barcode, print("Low"), print("Unknown"))
Short_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_child@meta.data) %in% Short_read_uninfected_child_medium_barcode & Short_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_uninfected_child@meta.data$Infection_tier)
Short_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_child@meta.data) %in% Short_read_uninfected_child_high_barcode & Short_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("High"),Short_read_uninfected_child@meta.data$Infection_tier)
Short_read_uninfected_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_uninfected_child@meta.data) %in% Short_read_uninfected_child_very_high_barcode & Short_read_uninfected_child@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_uninfected_child@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_uninfected_child@meta.data[Short_read_uninfected_child@meta.data$Infection_tier == "Unknown",]))


Short_read_uninfected_child <-importDemux(
    Short_read_uninfected_child,
    demuxlet.best = demuxlet_Short_read_uninfected_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Short_read_uninfected_child)
demux.calls.summary(Short_read_uninfected_child)



Short_read_VIC01_48hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_VIC01_48hpi_child/outs/filtered_feature_bc_matrix")
Short_read_VIC01_48hpi_child <- CreateSeuratObject(counts = Short_read_VIC01_48hpi_child)
demuxlet_Short_read_VIC01_48hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_VIC01_48hpi_child/outs/demuxlet/Short_read_VIC01_48hpi_child_nowarning.best", header=TRUE)

Short_read_VIC01_48hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_VIC01_48hpi_child/raw_feature_bc_matrix")
Short_read_VIC01_48hpi_child_v <- CreateSeuratObject(counts = Short_read_VIC01_48hpi_child_v)


#Select infected cells
Short_read_VIC01_48hpi_child_v@meta.data <- subset(Short_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_VIC01_48hpi_child_v_barcode <- rownames(Short_read_VIC01_48hpi_child_v@meta.data)
Short_read_VIC01_48hpi_child_v_barcode <- as.list(Short_read_VIC01_48hpi_child_v_barcode)

#Add infection level
Short_read_VIC01_48hpi_child@meta.data$Infection <- ifelse(rownames(Short_read_VIC01_48hpi_child@meta.data) %in% Short_read_VIC01_48hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection == "Uninfected",]))

Short_read_VIC01_48hpi_child_low_meta.data <- subset(Short_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_VIC01_48hpi_child_medium_meta.data <- subset(Short_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_VIC01_48hpi_child_high_meta.data <- subset(Short_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_VIC01_48hpi_child_very_high_meta.data <- subset(Short_read_VIC01_48hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_VIC01_48hpi_child_low_barcode <- rownames(Short_read_VIC01_48hpi_child_low_meta.data)
Short_read_VIC01_48hpi_child_low_barcode <- as.list(Short_read_VIC01_48hpi_child_low_barcode)

Short_read_VIC01_48hpi_child_medium_barcode <- rownames(Short_read_VIC01_48hpi_child_medium_meta.data)
Short_read_VIC01_48hpi_child_medium_barcode <- as.list(Short_read_VIC01_48hpi_child_medium_barcode)

Short_read_VIC01_48hpi_child_high_barcode <- rownames(Short_read_VIC01_48hpi_child_high_meta.data)
Short_read_VIC01_48hpi_child_high_barcode <- as.list(Short_read_VIC01_48hpi_child_high_barcode)

Short_read_VIC01_48hpi_child_very_high_barcode <- rownames(Short_read_VIC01_48hpi_child_very_high_meta.data)
Short_read_VIC01_48hpi_child_very_high_barcode <- as.list(Short_read_VIC01_48hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_child@meta.data) %in% Short_read_VIC01_48hpi_child_low_barcode, print("Low"), print("Unknown"))
Short_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_child@meta.data) %in% Short_read_VIC01_48hpi_child_medium_barcode & Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_VIC01_48hpi_child@meta.data$Infection_tier)
Short_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_child@meta.data) %in% Short_read_VIC01_48hpi_child_high_barcode & Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Short_read_VIC01_48hpi_child@meta.data$Infection_tier)
Short_read_VIC01_48hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_48hpi_child@meta.data) %in% Short_read_VIC01_48hpi_child_very_high_barcode & Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_VIC01_48hpi_child@meta.data$Infection_tier)



#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_VIC01_48hpi_child@meta.data[Short_read_VIC01_48hpi_child@meta.data$Infection_tier == "Unknown",]))


Short_read_VIC01_48hpi_child <-importDemux(
    Short_read_VIC01_48hpi_child,
    demuxlet.best = demuxlet_Short_read_VIC01_48hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)

demux.SNP.summary(Short_read_VIC01_48hpi_child)
demux.calls.summary(Short_read_VIC01_48hpi_child)


Short_read_VIC01_72hpi_child <- Read10X(data.dir = "/analysis/cellranger_new_run/run/Short_read_VIC01_72hpi_child/outs/filtered_feature_bc_matrix")
Short_read_VIC01_72hpi_child <- CreateSeuratObject(counts = Short_read_VIC01_72hpi_child)
demuxlet_Short_read_VIC01_72hpi_child <- read.table(file= "/analysis/cellranger_new_run/run/Short_read_VIC01_72hpi_child/outs/demuxlet/Short_read_VIC01_72hpi_child_nowarning.best", header=TRUE)

Short_read_VIC01_72hpi_child_v <- Read10X(data.dir = "/analysis/viral_barcode_matrix/Short_read_VIC01_72hpi_child/raw_feature_bc_matrix")
Short_read_VIC01_72hpi_child_v <- CreateSeuratObject(counts = Short_read_VIC01_72hpi_child_v)


#Select infected cells
Short_read_VIC01_72hpi_child_v@meta.data <- subset(Short_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10,)

#Isolate cell barcodes
Short_read_VIC01_72hpi_child_v_barcode <- rownames(Short_read_VIC01_72hpi_child_v@meta.data)
Short_read_VIC01_72hpi_child_v_barcode <- as.list(Short_read_VIC01_72hpi_child_v_barcode)

#Add infection level
Short_read_VIC01_72hpi_child@meta.data$Infection <- ifelse(rownames(Short_read_VIC01_72hpi_child@meta.data) %in% Short_read_VIC01_72hpi_child_v_barcode, print("Infected"),print("Uninfected"))


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection == "Infected",]))
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection == "Uninfected",]))

Short_read_VIC01_72hpi_child_low_meta.data <- subset(Short_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10 & nCount_RNA <100, )
Short_read_VIC01_72hpi_child_medium_meta.data <- subset(Short_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 100 & nCount_RNA <1000, )
Short_read_VIC01_72hpi_child_high_meta.data <- subset(Short_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 1000 & nCount_RNA <10000, )
Short_read_VIC01_72hpi_child_very_high_meta.data <- subset(Short_read_VIC01_72hpi_child_v@meta.data, nCount_RNA >= 10000, )


#Isolate cell barcodes
Short_read_VIC01_72hpi_child_low_barcode <- rownames(Short_read_VIC01_72hpi_child_low_meta.data)
Short_read_VIC01_72hpi_child_low_barcode <- as.list(Short_read_VIC01_72hpi_child_low_barcode)

Short_read_VIC01_72hpi_child_medium_barcode <- rownames(Short_read_VIC01_72hpi_child_medium_meta.data)
Short_read_VIC01_72hpi_child_medium_barcode <- as.list(Short_read_VIC01_72hpi_child_medium_barcode)

Short_read_VIC01_72hpi_child_high_barcode <- rownames(Short_read_VIC01_72hpi_child_high_meta.data)
Short_read_VIC01_72hpi_child_high_barcode <- as.list(Short_read_VIC01_72hpi_child_high_barcode)

Short_read_VIC01_72hpi_child_very_high_barcode <- rownames(Short_read_VIC01_72hpi_child_very_high_meta.data)
Short_read_VIC01_72hpi_child_very_high_barcode <- as.list(Short_read_VIC01_72hpi_child_very_high_barcode)
#-----------------------Final
#Add infection level
Short_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_child@meta.data) %in% Short_read_VIC01_72hpi_child_low_barcode, print("Low"), print("Unknown"))
Short_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_child@meta.data) %in% Short_read_VIC01_72hpi_child_medium_barcode & Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("Medium"),Short_read_VIC01_72hpi_child@meta.data$Infection_tier)
Short_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_child@meta.data) %in% Short_read_VIC01_72hpi_child_high_barcode & Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("High"),Short_read_VIC01_72hpi_child@meta.data$Infection_tier)
Short_read_VIC01_72hpi_child@meta.data$Infection_tier <- ifelse (rownames(Short_read_VIC01_72hpi_child@meta.data) %in% Short_read_VIC01_72hpi_child_very_high_barcode & Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown", print("Very High"),Short_read_VIC01_72hpi_child@meta.data$Infection_tier)


#Count number of cells with infected vs uninfected
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Low",]))
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Medium",]))
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "High",]))
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Very High",]))
print(nrow(Short_read_VIC01_72hpi_child@meta.data[Short_read_VIC01_72hpi_child@meta.data$Infection_tier == "Unknown",]))


Short_read_VIC01_72hpi_child <-importDemux(
    Short_read_VIC01_72hpi_child,
    demuxlet.best = demuxlet_Short_read_VIC01_72hpi_child,
    trim.before_ = TRUE,
    bypass.check = FALSE,
    verbose = TRUE
)


demux.SNP.summary(Short_read_VIC01_72hpi_child)
demux.calls.summary(Short_read_VIC01_72hpi_child)



# add metadata
Long_read_UK_72hpi_adult$type = "Long_read_UK_72hpi_adult"
Long_read_uninfected_adult$type = "Long_read_uninfected_adult"
Long_read_VIC01_48hpi_adult$type = "Long_read_VIC01_48hpi_adult"
Long_read_VIC01_72hpi_adult$type = "Long_read_VIC01_72hpi_adult"
Short_read_UK_72hpi_adult$type = "Short_read_UK_72hpi_adult"
Short_read_uninfected_adult$type = "Short_read_uninfected_adult"
Short_read_VIC01_48hpi_adult$type = "Short_read_VIC01_48hpi_adult"
Short_read_VIC01_72hpi_adult$type = "Short_read_VIC01_72hpi_adult"



Long_read_UK_72hpi_adult$group = "UK_72hpi_adult"
Long_read_uninfected_adult$group = "uninfected_adult"
Long_read_VIC01_48hpi_adult$group = "VIC01_48hpi_adult"
Long_read_VIC01_72hpi_adult$group = "VIC01_72hpi_adult"
Short_read_UK_72hpi_adult$group = "UK_72hpi_adult"
Short_read_uninfected_adult$group = "uninfected_adult"
Short_read_VIC01_48hpi_adult$group = "VIC01_48hpi_adult"
Short_read_VIC01_72hpi_adult$group = "VIC01_72hpi_adult"



Long_read_UK_72hpi_child$type = "Long_read_UK_72hpi_child"
Long_read_uninfected_child$type = "Long_read_uninfected_child"
Long_read_VIC01_48hpi_child$type = "Long_read_VIC01_48hpi_child"
Long_read_VIC01_72hpi_child$type = "Long_read_VIC01_72hpi_child"
Short_read_UK_72hpi_child$type = "Short_read_UK_72hpi_child"
Short_read_uninfected_child$type = "Short_read_uninfected_child"
Short_read_VIC01_48hpi_child$type = "Short_read_VIC01_48hpi_child"
Short_read_VIC01_72hpi_child$type = "Short_read_VIC01_72hpi_child"



Long_read_UK_72hpi_child$group = "UK_72hpi_child"
Long_read_uninfected_child$group = "uninfected_child"
Long_read_VIC01_48hpi_child$group = "VIC01_48hpi_child"
Long_read_VIC01_72hpi_child$group = "VIC01_72hpi_child"
Short_read_UK_72hpi_child$group = "UK_72hpi_child"
Short_read_uninfected_child$group = "uninfected_child"
Short_read_VIC01_48hpi_child$group = "VIC01_48hpi_child"
Short_read_VIC01_72hpi_child$group = "VIC01_72hpi_child"

#----------


#----------


# store mitochondrial percentage in object meta data
Long_read_UK_72hpi_adult <- PercentageFeatureSet(Long_read_UK_72hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Long_read_UK_72hpi_adult_f <- PercentageFeatureSet(Long_read_UK_72hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_UK_72hpi_adult <- PercentageFeatureSet(Long_read_UK_72hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_UK_72hpi_adult <- PercentageFeatureSet(Long_read_UK_72hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_UK_72hpi_adult <- SCTransform(Long_read_UK_72hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_UK_72hpi_adult_f <- SCTransform(Long_read_UK_72hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_UK_72hpi_adult_data <- subset(Long_read_UK_72hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_UK_72hpi_adult_selected_c <- WhichCells(Long_read_UK_72hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_UK_72hpi_adult_selected_f <- rownames(Long_read_UK_72hpi_adult_data)[Matrix::rowSums(Long_read_UK_72hpi_adult_data) > 3]

Long_read_UK_72hpi_adult_data.filt <- subset(Long_read_UK_72hpi_adult_data, features = Long_read_UK_72hpi_adult_selected_f, cells = Long_read_UK_72hpi_adult_selected_c)
dim(Long_read_UK_72hpi_adult_data.filt)




#mito/ribo filtering
Long_read_UK_72hpi_adult_selected_mito <- WhichCells(Long_read_UK_72hpi_adult_data.filt, expression = percent_mito < 20)
Long_read_UK_72hpi_adult_selected_ribo <- WhichCells(Long_read_UK_72hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_UK_72hpi_adult_data.filt <- subset(Long_read_UK_72hpi_adult_data.filt, cells = Long_read_UK_72hpi_adult_selected_mito)
Long_read_UK_72hpi_adult_data.filt <- subset(Long_read_UK_72hpi_adult_data.filt, cells = Long_read_UK_72hpi_adult_selected_ribo)

dim(Long_read_UK_72hpi_adult_data.filt)

table(Long_read_UK_72hpi_adult_data.filt$orig.ident)


Long_read_UK_72hpi_adult_data.filt_1 <- subset(Long_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Long_read_UK_72hpi_adult_data.filt_2 <- subset(Long_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Long_read_UK_72hpi_adult_data.filt_3 <- subset(Long_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_UK_72hpi_adult_data.filt_1 <- NormalizeData(Long_read_UK_72hpi_adult_data.filt_1)
Long_read_UK_72hpi_adult_data.filt_2 <- NormalizeData(Long_read_UK_72hpi_adult_data.filt_2)
Long_read_UK_72hpi_adult_data.filt_3 <- NormalizeData(Long_read_UK_72hpi_adult_data.filt_3)



Long_read_UK_72hpi_adult_data.filt_1 <- FindVariableFeatures(Long_read_UK_72hpi_adult_data.filt_1, selection.method = "vst")
Long_read_UK_72hpi_adult_data.filt_2 <- FindVariableFeatures(Long_read_UK_72hpi_adult_data.filt_2, selection.method = "vst")
Long_read_UK_72hpi_adult_data.filt_3 <- FindVariableFeatures(Long_read_UK_72hpi_adult_data.filt_3, selection.method = "vst")


Long_read_UK_72hpi_adult_data.filt_1 <- ScaleData(Long_read_UK_72hpi_adult_data.filt_1, features = rownames(Long_read_UK_72hpi_adult_data.filt_1))
Long_read_UK_72hpi_adult_data.filt_2 <- ScaleData(Long_read_UK_72hpi_adult_data.filt_2, features = rownames(Long_read_UK_72hpi_adult_data.filt_2))
Long_read_UK_72hpi_adult_data.filt_3 <- ScaleData(Long_read_UK_72hpi_adult_data.filt_3, features = rownames(Long_read_UK_72hpi_adult_data.filt_3))



Long_read_UK_72hpi_adult_data.filt_1 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_1, features = VariableFeatures(Long_read_UK_72hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_UK_72hpi_adult_data.filt_2 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_2, features = VariableFeatures(Long_read_UK_72hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_UK_72hpi_adult_data.filt_3 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_3, features = VariableFeatures(Long_read_UK_72hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)



Long_read_UK_72hpi_adult_data.filt_1 <- CellCycleScoring(Long_read_UK_72hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_UK_72hpi_adult_data.filt_2 <- CellCycleScoring(Long_read_UK_72hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_UK_72hpi_adult_data.filt_3 <- CellCycleScoring(Long_read_UK_72hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



#Alternate Workflow regression
Long_read_UK_72hpi_adult_data.filt_1$CC.Difference <- Long_read_UK_72hpi_adult_data.filt_1$S.Score - Long_read_UK_72hpi_adult_data.filt_1$G2M.Score
Long_read_UK_72hpi_adult_data.filt_2$CC.Difference <- Long_read_UK_72hpi_adult_data.filt_2$S.Score - Long_read_UK_72hpi_adult_data.filt_2$G2M.Score
Long_read_UK_72hpi_adult_data.filt_3$CC.Difference <- Long_read_UK_72hpi_adult_data.filt_3$S.Score - Long_read_UK_72hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_UK_72hpi_adult_data.filt_1 <- SCTransform(Long_read_UK_72hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_UK_72hpi_adult_data.filt_2 <- SCTransform(Long_read_UK_72hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_UK_72hpi_adult_data.filt_3 <- SCTransform(Long_read_UK_72hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)


Long_read_UK_72hpi_adult_data.filt_1 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_UK_72hpi_adult_data.filt_2 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_UK_72hpi_adult_data.filt_3 <- RunPCA(Long_read_UK_72hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Long_read_UK_72hpi_adult_data.filt_1)
DimPlot(Long_read_UK_72hpi_adult_data.filt_2)
DimPlot(Long_read_UK_72hpi_adult_data.filt_3)




# store mitochondrial percentage in object meta data
Long_read_uninfected_adult <- PercentageFeatureSet(Long_read_uninfected_adult, pattern = "^MT-", col.name = "percent_mito")
#Long_read_uninfected_adult_f <- PercentageFeatureSet(Long_read_uninfected_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_uninfected_adult <- PercentageFeatureSet(Long_read_uninfected_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_uninfected_adult <- PercentageFeatureSet(Long_read_uninfected_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_uninfected_adult <- SCTransform(Long_read_uninfected_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_uninfected_adult_f <- SCTransform(Long_read_uninfected_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_uninfected_adult_data <- subset(Long_read_uninfected_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_uninfected_adult_selected_c <- WhichCells(Long_read_uninfected_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_uninfected_adult_selected_f <- rownames(Long_read_uninfected_adult_data)[Matrix::rowSums(Long_read_uninfected_adult_data) > 3]

Long_read_uninfected_adult_data.filt <- subset(Long_read_uninfected_adult_data, features = Long_read_uninfected_adult_selected_f, cells = Long_read_uninfected_adult_selected_c)
dim(Long_read_uninfected_adult_data.filt)




#mito/ribo filtering
Long_read_uninfected_adult_selected_mito <- WhichCells(Long_read_uninfected_adult_data.filt, expression = percent_mito < 20)
Long_read_uninfected_adult_selected_ribo <- WhichCells(Long_read_uninfected_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_uninfected_adult_data.filt <- subset(Long_read_uninfected_adult_data.filt, cells = Long_read_uninfected_adult_selected_mito)
Long_read_uninfected_adult_data.filt <- subset(Long_read_uninfected_adult_data.filt, cells = Long_read_uninfected_adult_selected_ribo)

dim(Long_read_uninfected_adult_data.filt)

table(Long_read_uninfected_adult_data.filt$orig.ident)


Long_read_uninfected_adult_data.filt_1 <- subset(Long_read_uninfected_adult_data.filt, subset = Sample =="Sample1_Sample1")
Long_read_uninfected_adult_data.filt_2 <- subset(Long_read_uninfected_adult_data.filt, subset = Sample =="Sample2_Sample2")
Long_read_uninfected_adult_data.filt_3 <- subset(Long_read_uninfected_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_uninfected_adult_data.filt_1 <- NormalizeData(Long_read_uninfected_adult_data.filt_1)
Long_read_uninfected_adult_data.filt_2 <- NormalizeData(Long_read_uninfected_adult_data.filt_2)
Long_read_uninfected_adult_data.filt_3 <- NormalizeData(Long_read_uninfected_adult_data.filt_3)



Long_read_uninfected_adult_data.filt_1 <- FindVariableFeatures(Long_read_uninfected_adult_data.filt_1, selection.method = "vst")
Long_read_uninfected_adult_data.filt_2 <- FindVariableFeatures(Long_read_uninfected_adult_data.filt_2, selection.method = "vst")
Long_read_uninfected_adult_data.filt_3 <- FindVariableFeatures(Long_read_uninfected_adult_data.filt_3, selection.method = "vst")


Long_read_uninfected_adult_data.filt_1 <- ScaleData(Long_read_uninfected_adult_data.filt_1, features = rownames(Long_read_uninfected_adult_data.filt_1))
Long_read_uninfected_adult_data.filt_2 <- ScaleData(Long_read_uninfected_adult_data.filt_2, features = rownames(Long_read_uninfected_adult_data.filt_2))
Long_read_uninfected_adult_data.filt_3 <- ScaleData(Long_read_uninfected_adult_data.filt_3, features = rownames(Long_read_uninfected_adult_data.filt_3))



Long_read_uninfected_adult_data.filt_1 <- RunPCA(Long_read_uninfected_adult_data.filt_1, features = VariableFeatures(Long_read_uninfected_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_uninfected_adult_data.filt_2 <- RunPCA(Long_read_uninfected_adult_data.filt_2, features = VariableFeatures(Long_read_uninfected_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_uninfected_adult_data.filt_3 <- RunPCA(Long_read_uninfected_adult_data.filt_3, features = VariableFeatures(Long_read_uninfected_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)




Long_read_uninfected_adult_data.filt_1 <- CellCycleScoring(Long_read_uninfected_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_uninfected_adult_data.filt_2 <- CellCycleScoring(Long_read_uninfected_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_uninfected_adult_data.filt_3 <- CellCycleScoring(Long_read_uninfected_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



#Alternate Workflow regression
Long_read_uninfected_adult_data.filt_1$CC.Difference <- Long_read_uninfected_adult_data.filt_1$S.Score - Long_read_uninfected_adult_data.filt_1$G2M.Score
Long_read_uninfected_adult_data.filt_2$CC.Difference <- Long_read_uninfected_adult_data.filt_2$S.Score - Long_read_uninfected_adult_data.filt_2$G2M.Score
Long_read_uninfected_adult_data.filt_3$CC.Difference <- Long_read_uninfected_adult_data.filt_3$S.Score - Long_read_uninfected_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_uninfected_adult_data.filt_1 <- SCTransform(Long_read_uninfected_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_uninfected_adult_data.filt_2 <- SCTransform(Long_read_uninfected_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_uninfected_adult_data.filt_3 <- SCTransform(Long_read_uninfected_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)


Long_read_uninfected_adult_data.filt_1 <- RunPCA(Long_read_uninfected_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_uninfected_adult_data.filt_2 <- RunPCA(Long_read_uninfected_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_uninfected_adult_data.filt_3 <- RunPCA(Long_read_uninfected_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Long_read_uninfected_adult_data.filt_1)
DimPlot(Long_read_uninfected_adult_data.filt_2)
DimPlot(Long_read_uninfected_adult_data.filt_3)




# store mitochondrial percentage in object meta data
Long_read_VIC01_48hpi_adult <- PercentageFeatureSet(Long_read_VIC01_48hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Long_read_VIC01_48hpi_adult_f <- PercentageFeatureSet(Long_read_VIC01_48hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_VIC01_48hpi_adult <- PercentageFeatureSet(Long_read_VIC01_48hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_VIC01_48hpi_adult <- PercentageFeatureSet(Long_read_VIC01_48hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_VIC01_48hpi_adult <- SCTransform(Long_read_VIC01_48hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_VIC01_48hpi_adult_f <- SCTransform(Long_read_VIC01_48hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_VIC01_48hpi_adult_data <- subset(Long_read_VIC01_48hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_VIC01_48hpi_adult_selected_c <- WhichCells(Long_read_VIC01_48hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_VIC01_48hpi_adult_selected_f <- rownames(Long_read_VIC01_48hpi_adult_data)[Matrix::rowSums(Long_read_VIC01_48hpi_adult_data) > 3]

Long_read_VIC01_48hpi_adult_data.filt <- subset(Long_read_VIC01_48hpi_adult_data, features = Long_read_VIC01_48hpi_adult_selected_f, cells = Long_read_VIC01_48hpi_adult_selected_c)
dim(Long_read_VIC01_48hpi_adult_data.filt)




#mito/ribo filtering
Long_read_VIC01_48hpi_adult_selected_mito <- WhichCells(Long_read_VIC01_48hpi_adult_data.filt, expression = percent_mito < 20)
Long_read_VIC01_48hpi_adult_selected_ribo <- WhichCells(Long_read_VIC01_48hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_VIC01_48hpi_adult_data.filt <- subset(Long_read_VIC01_48hpi_adult_data.filt, cells = Long_read_VIC01_48hpi_adult_selected_mito)
Long_read_VIC01_48hpi_adult_data.filt <- subset(Long_read_VIC01_48hpi_adult_data.filt, cells = Long_read_VIC01_48hpi_adult_selected_ribo)

dim(Long_read_VIC01_48hpi_adult_data.filt)

table(Long_read_VIC01_48hpi_adult_data.filt$orig.ident)


Long_read_VIC01_48hpi_adult_data.filt_1 <- subset(Long_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Long_read_VIC01_48hpi_adult_data.filt_2 <- subset(Long_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Long_read_VIC01_48hpi_adult_data.filt_3 <- subset(Long_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_VIC01_48hpi_adult_data.filt_1 <- NormalizeData(Long_read_VIC01_48hpi_adult_data.filt_1)
Long_read_VIC01_48hpi_adult_data.filt_2 <- NormalizeData(Long_read_VIC01_48hpi_adult_data.filt_2)
Long_read_VIC01_48hpi_adult_data.filt_3 <- NormalizeData(Long_read_VIC01_48hpi_adult_data.filt_3)



Long_read_VIC01_48hpi_adult_data.filt_1 <- FindVariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_1, selection.method = "vst")
Long_read_VIC01_48hpi_adult_data.filt_2 <- FindVariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_2, selection.method = "vst")
Long_read_VIC01_48hpi_adult_data.filt_3 <- FindVariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_3, selection.method = "vst")


Long_read_VIC01_48hpi_adult_data.filt_1 <- ScaleData(Long_read_VIC01_48hpi_adult_data.filt_1, features = rownames(Long_read_VIC01_48hpi_adult_data.filt_1))
Long_read_VIC01_48hpi_adult_data.filt_2 <- ScaleData(Long_read_VIC01_48hpi_adult_data.filt_2, features = rownames(Long_read_VIC01_48hpi_adult_data.filt_2))
Long_read_VIC01_48hpi_adult_data.filt_3 <- ScaleData(Long_read_VIC01_48hpi_adult_data.filt_3, features = rownames(Long_read_VIC01_48hpi_adult_data.filt_3))



Long_read_VIC01_48hpi_adult_data.filt_1 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_1, features = VariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_48hpi_adult_data.filt_2 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_2, features = VariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_48hpi_adult_data.filt_3 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_3, features = VariableFeatures(Long_read_VIC01_48hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)




Long_read_VIC01_48hpi_adult_data.filt_1 <- CellCycleScoring(Long_read_VIC01_48hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_48hpi_adult_data.filt_2 <- CellCycleScoring(Long_read_VIC01_48hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_48hpi_adult_data.filt_3 <- CellCycleScoring(Long_read_VIC01_48hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Long_read_VIC01_48hpi_adult_data.filt_1$CC.Difference <- Long_read_VIC01_48hpi_adult_data.filt_1$S.Score - Long_read_VIC01_48hpi_adult_data.filt_1$G2M.Score
Long_read_VIC01_48hpi_adult_data.filt_2$CC.Difference <- Long_read_VIC01_48hpi_adult_data.filt_2$S.Score - Long_read_VIC01_48hpi_adult_data.filt_2$G2M.Score
Long_read_VIC01_48hpi_adult_data.filt_3$CC.Difference <- Long_read_VIC01_48hpi_adult_data.filt_3$S.Score - Long_read_VIC01_48hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_VIC01_48hpi_adult_data.filt_1 <- SCTransform(Long_read_VIC01_48hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_48hpi_adult_data.filt_2 <- SCTransform(Long_read_VIC01_48hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_48hpi_adult_data.filt_3 <- SCTransform(Long_read_VIC01_48hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_VIC01_48hpi_adult_data.filt_1 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_48hpi_adult_data.filt_2 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_48hpi_adult_data.filt_3 <- RunPCA(Long_read_VIC01_48hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Long_read_VIC01_48hpi_adult_data.filt_1)
DimPlot(Long_read_VIC01_48hpi_adult_data.filt_2)
DimPlot(Long_read_VIC01_48hpi_adult_data.filt_3)



# store mitochondrial percentage in object meta data
Long_read_VIC01_72hpi_adult <- PercentageFeatureSet(Long_read_VIC01_72hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Long_read_VIC01_72hpi_adult_f <- PercentageFeatureSet(Long_read_VIC01_72hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_VIC01_72hpi_adult <- PercentageFeatureSet(Long_read_VIC01_72hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_VIC01_72hpi_adult <- PercentageFeatureSet(Long_read_VIC01_72hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_VIC01_72hpi_adult <- SCTransform(Long_read_VIC01_72hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_VIC01_72hpi_adult_f <- SCTransform(Long_read_VIC01_72hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_VIC01_72hpi_adult_data <- subset(Long_read_VIC01_72hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_VIC01_72hpi_adult_selected_c <- WhichCells(Long_read_VIC01_72hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_VIC01_72hpi_adult_selected_f <- rownames(Long_read_VIC01_72hpi_adult_data)[Matrix::rowSums(Long_read_VIC01_72hpi_adult_data) > 3]

Long_read_VIC01_72hpi_adult_data.filt <- subset(Long_read_VIC01_72hpi_adult_data, features = Long_read_VIC01_72hpi_adult_selected_f, cells = Long_read_VIC01_72hpi_adult_selected_c)
dim(Long_read_VIC01_72hpi_adult_data.filt)




#mito/ribo filtering
Long_read_VIC01_72hpi_adult_selected_mito <- WhichCells(Long_read_VIC01_72hpi_adult_data.filt, expression = percent_mito < 20)
Long_read_VIC01_72hpi_adult_selected_ribo <- WhichCells(Long_read_VIC01_72hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_VIC01_72hpi_adult_data.filt <- subset(Long_read_VIC01_72hpi_adult_data.filt, cells = Long_read_VIC01_72hpi_adult_selected_mito)
Long_read_VIC01_72hpi_adult_data.filt <- subset(Long_read_VIC01_72hpi_adult_data.filt, cells = Long_read_VIC01_72hpi_adult_selected_ribo)

dim(Long_read_VIC01_72hpi_adult_data.filt)

table(Long_read_VIC01_72hpi_adult_data.filt$orig.ident)


Long_read_VIC01_72hpi_adult_data.filt_1 <- subset(Long_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Long_read_VIC01_72hpi_adult_data.filt_2 <- subset(Long_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Long_read_VIC01_72hpi_adult_data.filt_3 <- subset(Long_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_VIC01_72hpi_adult_data.filt_1 <- NormalizeData(Long_read_VIC01_72hpi_adult_data.filt_1)
Long_read_VIC01_72hpi_adult_data.filt_2 <- NormalizeData(Long_read_VIC01_72hpi_adult_data.filt_2)
Long_read_VIC01_72hpi_adult_data.filt_3 <- NormalizeData(Long_read_VIC01_72hpi_adult_data.filt_3)



Long_read_VIC01_72hpi_adult_data.filt_1 <- FindVariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_1, selection.method = "vst")
Long_read_VIC01_72hpi_adult_data.filt_2 <- FindVariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_2, selection.method = "vst")
Long_read_VIC01_72hpi_adult_data.filt_3 <- FindVariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_3, selection.method = "vst")


Long_read_VIC01_72hpi_adult_data.filt_1 <- ScaleData(Long_read_VIC01_72hpi_adult_data.filt_1, features = rownames(Long_read_VIC01_72hpi_adult_data.filt_1))
Long_read_VIC01_72hpi_adult_data.filt_2 <- ScaleData(Long_read_VIC01_72hpi_adult_data.filt_2, features = rownames(Long_read_VIC01_72hpi_adult_data.filt_2))
Long_read_VIC01_72hpi_adult_data.filt_3 <- ScaleData(Long_read_VIC01_72hpi_adult_data.filt_3, features = rownames(Long_read_VIC01_72hpi_adult_data.filt_3))



Long_read_VIC01_72hpi_adult_data.filt_1 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_1, features = VariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_72hpi_adult_data.filt_2 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_2, features = VariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_72hpi_adult_data.filt_3 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_3, features = VariableFeatures(Long_read_VIC01_72hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)



Long_read_VIC01_72hpi_adult_data.filt_1 <- CellCycleScoring(Long_read_VIC01_72hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_72hpi_adult_data.filt_2 <- CellCycleScoring(Long_read_VIC01_72hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_72hpi_adult_data.filt_3 <- CellCycleScoring(Long_read_VIC01_72hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Long_read_VIC01_72hpi_adult_data.filt_1$CC.Difference <- Long_read_VIC01_72hpi_adult_data.filt_1$S.Score - Long_read_VIC01_72hpi_adult_data.filt_1$G2M.Score
Long_read_VIC01_72hpi_adult_data.filt_2$CC.Difference <- Long_read_VIC01_72hpi_adult_data.filt_2$S.Score - Long_read_VIC01_72hpi_adult_data.filt_2$G2M.Score
Long_read_VIC01_72hpi_adult_data.filt_3$CC.Difference <- Long_read_VIC01_72hpi_adult_data.filt_3$S.Score - Long_read_VIC01_72hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_VIC01_72hpi_adult_data.filt_1 <- SCTransform(Long_read_VIC01_72hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_72hpi_adult_data.filt_2 <- SCTransform(Long_read_VIC01_72hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_72hpi_adult_data.filt_3 <- SCTransform(Long_read_VIC01_72hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_VIC01_72hpi_adult_data.filt_1 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_72hpi_adult_data.filt_2 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_72hpi_adult_data.filt_3 <- RunPCA(Long_read_VIC01_72hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)



DimPlot(Long_read_VIC01_72hpi_adult_data.filt_1)
DimPlot(Long_read_VIC01_72hpi_adult_data.filt_2)
DimPlot(Long_read_VIC01_72hpi_adult_data.filt_3)



# store mitochondrial percentage in object meta data
Short_read_UK_72hpi_adult <- PercentageFeatureSet(Short_read_UK_72hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Short_read_UK_72hpi_adult_f <- PercentageFeatureSet(Short_read_UK_72hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_UK_72hpi_adult <- PercentageFeatureSet(Short_read_UK_72hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_UK_72hpi_adult <- PercentageFeatureSet(Short_read_UK_72hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_UK_72hpi_adult <- SCTransform(Short_read_UK_72hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_UK_72hpi_adult_f <- SCTransform(Short_read_UK_72hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_UK_72hpi_adult_data <- subset(Short_read_UK_72hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_UK_72hpi_adult_selected_c <- WhichCells(Short_read_UK_72hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_UK_72hpi_adult_selected_f <- rownames(Short_read_UK_72hpi_adult_data)[Matrix::rowSums(Short_read_UK_72hpi_adult_data) > 3]

Short_read_UK_72hpi_adult_data.filt <- subset(Short_read_UK_72hpi_adult_data, features = Short_read_UK_72hpi_adult_selected_f, cells = Short_read_UK_72hpi_adult_selected_c)
dim(Short_read_UK_72hpi_adult_data.filt)




#mito/ribo filtering
Short_read_UK_72hpi_adult_selected_mito <- WhichCells(Short_read_UK_72hpi_adult_data.filt, expression = percent_mito < 20)
Short_read_UK_72hpi_adult_selected_ribo <- WhichCells(Short_read_UK_72hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_UK_72hpi_adult_data.filt <- subset(Short_read_UK_72hpi_adult_data.filt, cells = Short_read_UK_72hpi_adult_selected_mito)
Short_read_UK_72hpi_adult_data.filt <- subset(Short_read_UK_72hpi_adult_data.filt, cells = Short_read_UK_72hpi_adult_selected_ribo)

dim(Short_read_UK_72hpi_adult_data.filt)

table(Short_read_UK_72hpi_adult_data.filt$orig.ident)


Short_read_UK_72hpi_adult_data.filt_1 <- subset(Short_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Short_read_UK_72hpi_adult_data.filt_2 <- subset(Short_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Short_read_UK_72hpi_adult_data.filt_3 <- subset(Short_read_UK_72hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_UK_72hpi_adult_data.filt_1 <- NormalizeData(Short_read_UK_72hpi_adult_data.filt_1)
Short_read_UK_72hpi_adult_data.filt_2 <- NormalizeData(Short_read_UK_72hpi_adult_data.filt_2)
Short_read_UK_72hpi_adult_data.filt_3 <- NormalizeData(Short_read_UK_72hpi_adult_data.filt_3)



Short_read_UK_72hpi_adult_data.filt_1 <- FindVariableFeatures(Short_read_UK_72hpi_adult_data.filt_1, selection.method = "vst")
Short_read_UK_72hpi_adult_data.filt_2 <- FindVariableFeatures(Short_read_UK_72hpi_adult_data.filt_2, selection.method = "vst")
Short_read_UK_72hpi_adult_data.filt_3 <- FindVariableFeatures(Short_read_UK_72hpi_adult_data.filt_3, selection.method = "vst")


Short_read_UK_72hpi_adult_data.filt_1 <- ScaleData(Short_read_UK_72hpi_adult_data.filt_1, features = rownames(Short_read_UK_72hpi_adult_data.filt_1))
Short_read_UK_72hpi_adult_data.filt_2 <- ScaleData(Short_read_UK_72hpi_adult_data.filt_2, features = rownames(Short_read_UK_72hpi_adult_data.filt_2))
Short_read_UK_72hpi_adult_data.filt_3 <- ScaleData(Short_read_UK_72hpi_adult_data.filt_3, features = rownames(Short_read_UK_72hpi_adult_data.filt_3))



Short_read_UK_72hpi_adult_data.filt_1 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_1, features = VariableFeatures(Short_read_UK_72hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_UK_72hpi_adult_data.filt_2 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_2, features = VariableFeatures(Short_read_UK_72hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_UK_72hpi_adult_data.filt_3 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_3, features = VariableFeatures(Short_read_UK_72hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)


Short_read_UK_72hpi_adult_data.filt_1 <- CellCycleScoring(Short_read_UK_72hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_UK_72hpi_adult_data.filt_2 <- CellCycleScoring(Short_read_UK_72hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_UK_72hpi_adult_data.filt_3 <- CellCycleScoring(Short_read_UK_72hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#Alternate Workflow regression
Short_read_UK_72hpi_adult_data.filt_1$CC.Difference <- Short_read_UK_72hpi_adult_data.filt_1$S.Score - Short_read_UK_72hpi_adult_data.filt_1$G2M.Score
Short_read_UK_72hpi_adult_data.filt_2$CC.Difference <- Short_read_UK_72hpi_adult_data.filt_2$S.Score - Short_read_UK_72hpi_adult_data.filt_2$G2M.Score
Short_read_UK_72hpi_adult_data.filt_3$CC.Difference <- Short_read_UK_72hpi_adult_data.filt_3$S.Score - Short_read_UK_72hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_UK_72hpi_adult_data.filt_1 <- SCTransform(Short_read_UK_72hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_UK_72hpi_adult_data.filt_2 <- SCTransform(Short_read_UK_72hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_UK_72hpi_adult_data.filt_3 <- SCTransform(Short_read_UK_72hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_UK_72hpi_adult_data.filt_1 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_UK_72hpi_adult_data.filt_2 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_UK_72hpi_adult_data.filt_3 <- RunPCA(Short_read_UK_72hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_UK_72hpi_adult_data.filt_1)
DimPlot(Short_read_UK_72hpi_adult_data.filt_2)
DimPlot(Short_read_UK_72hpi_adult_data.filt_3)



# store mitochondrial percentage in object meta data
Short_read_uninfected_adult <- PercentageFeatureSet(Short_read_uninfected_adult, pattern = "^MT-", col.name = "percent_mito")
#Short_read_uninfected_adult_f <- PercentageFeatureSet(Short_read_uninfected_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_uninfected_adult <- PercentageFeatureSet(Short_read_uninfected_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_uninfected_adult <- PercentageFeatureSet(Short_read_uninfected_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_uninfected_adult <- SCTransform(Short_read_uninfected_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_uninfected_adult_f <- SCTransform(Short_read_uninfected_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_uninfected_adult_data <- subset(Short_read_uninfected_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_uninfected_adult_selected_c <- WhichCells(Short_read_uninfected_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_uninfected_adult_selected_f <- rownames(Short_read_uninfected_adult_data)[Matrix::rowSums(Short_read_uninfected_adult_data) > 3]

Short_read_uninfected_adult_data.filt <- subset(Short_read_uninfected_adult_data, features = Short_read_uninfected_adult_selected_f, cells = Short_read_uninfected_adult_selected_c)
dim(Short_read_uninfected_adult_data.filt)




#mito/ribo filtering
Short_read_uninfected_adult_selected_mito <- WhichCells(Short_read_uninfected_adult_data.filt, expression = percent_mito < 20)
Short_read_uninfected_adult_selected_ribo <- WhichCells(Short_read_uninfected_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_uninfected_adult_data.filt <- subset(Short_read_uninfected_adult_data.filt, cells = Short_read_uninfected_adult_selected_mito)
Short_read_uninfected_adult_data.filt <- subset(Short_read_uninfected_adult_data.filt, cells = Short_read_uninfected_adult_selected_ribo)

dim(Short_read_uninfected_adult_data.filt)

table(Short_read_uninfected_adult_data.filt$orig.ident)


Short_read_uninfected_adult_data.filt_1 <- subset(Short_read_uninfected_adult_data.filt, subset = Sample =="Sample1_Sample1")
Short_read_uninfected_adult_data.filt_2 <- subset(Short_read_uninfected_adult_data.filt, subset = Sample =="Sample2_Sample2")
Short_read_uninfected_adult_data.filt_3 <- subset(Short_read_uninfected_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_uninfected_adult_data.filt_1 <- NormalizeData(Short_read_uninfected_adult_data.filt_1)
Short_read_uninfected_adult_data.filt_2 <- NormalizeData(Short_read_uninfected_adult_data.filt_2)
Short_read_uninfected_adult_data.filt_3 <- NormalizeData(Short_read_uninfected_adult_data.filt_3)



Short_read_uninfected_adult_data.filt_1 <- FindVariableFeatures(Short_read_uninfected_adult_data.filt_1, selection.method = "vst")
Short_read_uninfected_adult_data.filt_2 <- FindVariableFeatures(Short_read_uninfected_adult_data.filt_2, selection.method = "vst")
Short_read_uninfected_adult_data.filt_3 <- FindVariableFeatures(Short_read_uninfected_adult_data.filt_3, selection.method = "vst")


Short_read_uninfected_adult_data.filt_1 <- ScaleData(Short_read_uninfected_adult_data.filt_1, features = rownames(Short_read_uninfected_adult_data.filt_1))
Short_read_uninfected_adult_data.filt_2 <- ScaleData(Short_read_uninfected_adult_data.filt_2, features = rownames(Short_read_uninfected_adult_data.filt_2))
Short_read_uninfected_adult_data.filt_3 <- ScaleData(Short_read_uninfected_adult_data.filt_3, features = rownames(Short_read_uninfected_adult_data.filt_3))



Short_read_uninfected_adult_data.filt_1 <- RunPCA(Short_read_uninfected_adult_data.filt_1, features = VariableFeatures(Short_read_uninfected_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_uninfected_adult_data.filt_2 <- RunPCA(Short_read_uninfected_adult_data.filt_2, features = VariableFeatures(Short_read_uninfected_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_uninfected_adult_data.filt_3 <- RunPCA(Short_read_uninfected_adult_data.filt_3, features = VariableFeatures(Short_read_uninfected_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)



Short_read_uninfected_adult_data.filt_1 <- CellCycleScoring(Short_read_uninfected_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_uninfected_adult_data.filt_2 <- CellCycleScoring(Short_read_uninfected_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_uninfected_adult_data.filt_3 <- CellCycleScoring(Short_read_uninfected_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)





#Alternate Workflow regression
Short_read_uninfected_adult_data.filt_1$CC.Difference <- Short_read_uninfected_adult_data.filt_1$S.Score - Short_read_uninfected_adult_data.filt_1$G2M.Score
Short_read_uninfected_adult_data.filt_2$CC.Difference <- Short_read_uninfected_adult_data.filt_2$S.Score - Short_read_uninfected_adult_data.filt_2$G2M.Score
Short_read_uninfected_adult_data.filt_3$CC.Difference <- Short_read_uninfected_adult_data.filt_3$S.Score - Short_read_uninfected_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_uninfected_adult_data.filt_1 <- SCTransform(Short_read_uninfected_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_uninfected_adult_data.filt_2 <- SCTransform(Short_read_uninfected_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_uninfected_adult_data.filt_3 <- SCTransform(Short_read_uninfected_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_uninfected_adult_data.filt_1 <- RunPCA(Short_read_uninfected_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_uninfected_adult_data.filt_2 <- RunPCA(Short_read_uninfected_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_uninfected_adult_data.filt_3 <- RunPCA(Short_read_uninfected_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_uninfected_adult_data.filt_1)
DimPlot(Short_read_uninfected_adult_data.filt_2)
DimPlot(Short_read_uninfected_adult_data.filt_3)


# store mitochondrial percentage in object meta data
Short_read_VIC01_48hpi_adult <- PercentageFeatureSet(Short_read_VIC01_48hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Short_read_VIC01_48hpi_adult_f <- PercentageFeatureSet(Short_read_VIC01_48hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_VIC01_48hpi_adult <- PercentageFeatureSet(Short_read_VIC01_48hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_VIC01_48hpi_adult <- PercentageFeatureSet(Short_read_VIC01_48hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_VIC01_48hpi_adult <- SCTransform(Short_read_VIC01_48hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_VIC01_48hpi_adult_f <- SCTransform(Short_read_VIC01_48hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_VIC01_48hpi_adult_data <- subset(Short_read_VIC01_48hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_VIC01_48hpi_adult_selected_c <- WhichCells(Short_read_VIC01_48hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_VIC01_48hpi_adult_selected_f <- rownames(Short_read_VIC01_48hpi_adult_data)[Matrix::rowSums(Short_read_VIC01_48hpi_adult_data) > 3]

Short_read_VIC01_48hpi_adult_data.filt <- subset(Short_read_VIC01_48hpi_adult_data, features = Short_read_VIC01_48hpi_adult_selected_f, cells = Short_read_VIC01_48hpi_adult_selected_c)
dim(Short_read_VIC01_48hpi_adult_data.filt)




#mito/ribo filtering
Short_read_VIC01_48hpi_adult_selected_mito <- WhichCells(Short_read_VIC01_48hpi_adult_data.filt, expression = percent_mito < 20)
Short_read_VIC01_48hpi_adult_selected_ribo <- WhichCells(Short_read_VIC01_48hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_VIC01_48hpi_adult_data.filt <- subset(Short_read_VIC01_48hpi_adult_data.filt, cells = Short_read_VIC01_48hpi_adult_selected_mito)
Short_read_VIC01_48hpi_adult_data.filt <- subset(Short_read_VIC01_48hpi_adult_data.filt, cells = Short_read_VIC01_48hpi_adult_selected_ribo)

dim(Short_read_VIC01_48hpi_adult_data.filt)

table(Short_read_VIC01_48hpi_adult_data.filt$orig.ident)


Short_read_VIC01_48hpi_adult_data.filt_1 <- subset(Short_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Short_read_VIC01_48hpi_adult_data.filt_2 <- subset(Short_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Short_read_VIC01_48hpi_adult_data.filt_3 <- subset(Short_read_VIC01_48hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_VIC01_48hpi_adult_data.filt_1 <- NormalizeData(Short_read_VIC01_48hpi_adult_data.filt_1)
Short_read_VIC01_48hpi_adult_data.filt_2 <- NormalizeData(Short_read_VIC01_48hpi_adult_data.filt_2)
Short_read_VIC01_48hpi_adult_data.filt_3 <- NormalizeData(Short_read_VIC01_48hpi_adult_data.filt_3)



Short_read_VIC01_48hpi_adult_data.filt_1 <- FindVariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_1, selection.method = "vst")
Short_read_VIC01_48hpi_adult_data.filt_2 <- FindVariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_2, selection.method = "vst")
Short_read_VIC01_48hpi_adult_data.filt_3 <- FindVariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_3, selection.method = "vst")


Short_read_VIC01_48hpi_adult_data.filt_1 <- ScaleData(Short_read_VIC01_48hpi_adult_data.filt_1, features = rownames(Short_read_VIC01_48hpi_adult_data.filt_1))
Short_read_VIC01_48hpi_adult_data.filt_2 <- ScaleData(Short_read_VIC01_48hpi_adult_data.filt_2, features = rownames(Short_read_VIC01_48hpi_adult_data.filt_2))
Short_read_VIC01_48hpi_adult_data.filt_3 <- ScaleData(Short_read_VIC01_48hpi_adult_data.filt_3, features = rownames(Short_read_VIC01_48hpi_adult_data.filt_3))



Short_read_VIC01_48hpi_adult_data.filt_1 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_1, features = VariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_48hpi_adult_data.filt_2 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_2, features = VariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_48hpi_adult_data.filt_3 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_3, features = VariableFeatures(Short_read_VIC01_48hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)


Short_read_VIC01_48hpi_adult_data.filt_1 <- CellCycleScoring(Short_read_VIC01_48hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_48hpi_adult_data.filt_2 <- CellCycleScoring(Short_read_VIC01_48hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_48hpi_adult_data.filt_3 <- CellCycleScoring(Short_read_VIC01_48hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#Alternate Workflow regression
Short_read_VIC01_48hpi_adult_data.filt_1$CC.Difference <- Short_read_VIC01_48hpi_adult_data.filt_1$S.Score - Short_read_VIC01_48hpi_adult_data.filt_1$G2M.Score
Short_read_VIC01_48hpi_adult_data.filt_2$CC.Difference <- Short_read_VIC01_48hpi_adult_data.filt_2$S.Score - Short_read_VIC01_48hpi_adult_data.filt_2$G2M.Score
Short_read_VIC01_48hpi_adult_data.filt_3$CC.Difference <- Short_read_VIC01_48hpi_adult_data.filt_3$S.Score - Short_read_VIC01_48hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_VIC01_48hpi_adult_data.filt_1 <- SCTransform(Short_read_VIC01_48hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_48hpi_adult_data.filt_2 <- SCTransform(Short_read_VIC01_48hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_48hpi_adult_data.filt_3 <- SCTransform(Short_read_VIC01_48hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_VIC01_48hpi_adult_data.filt_1 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_48hpi_adult_data.filt_2 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_48hpi_adult_data.filt_3 <- RunPCA(Short_read_VIC01_48hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_VIC01_48hpi_adult_data.filt_1)
DimPlot(Short_read_VIC01_48hpi_adult_data.filt_2)
DimPlot(Short_read_VIC01_48hpi_adult_data.filt_3)



# store mitochondrial percentage in object meta data
Short_read_VIC01_72hpi_adult <- PercentageFeatureSet(Short_read_VIC01_72hpi_adult, pattern = "^MT-", col.name = "percent_mito")
#Short_read_VIC01_72hpi_adult_f <- PercentageFeatureSet(Short_read_VIC01_72hpi_adult_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_VIC01_72hpi_adult <- PercentageFeatureSet(Short_read_VIC01_72hpi_adult, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_VIC01_72hpi_adult <- PercentageFeatureSet(Short_read_VIC01_72hpi_adult, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_VIC01_72hpi_adult <- SCTransform(Short_read_VIC01_72hpi_adult, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_VIC01_72hpi_adult_f <- SCTransform(Short_read_VIC01_72hpi_adult_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_VIC01_72hpi_adult_data <- subset(Short_read_VIC01_72hpi_adult, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_VIC01_72hpi_adult_selected_c <- WhichCells(Short_read_VIC01_72hpi_adult_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_VIC01_72hpi_adult_selected_f <- rownames(Short_read_VIC01_72hpi_adult_data)[Matrix::rowSums(Short_read_VIC01_72hpi_adult_data) > 3]

Short_read_VIC01_72hpi_adult_data.filt <- subset(Short_read_VIC01_72hpi_adult_data, features = Short_read_VIC01_72hpi_adult_selected_f, cells = Short_read_VIC01_72hpi_adult_selected_c)
dim(Short_read_VIC01_72hpi_adult_data.filt)




#mito/ribo filtering
Short_read_VIC01_72hpi_adult_selected_mito <- WhichCells(Short_read_VIC01_72hpi_adult_data.filt, expression = percent_mito < 20)
Short_read_VIC01_72hpi_adult_selected_ribo <- WhichCells(Short_read_VIC01_72hpi_adult_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_VIC01_72hpi_adult_data.filt <- subset(Short_read_VIC01_72hpi_adult_data.filt, cells = Short_read_VIC01_72hpi_adult_selected_mito)
Short_read_VIC01_72hpi_adult_data.filt <- subset(Short_read_VIC01_72hpi_adult_data.filt, cells = Short_read_VIC01_72hpi_adult_selected_ribo)

dim(Short_read_VIC01_72hpi_adult_data.filt)

table(Short_read_VIC01_72hpi_adult_data.filt$orig.ident)


Short_read_VIC01_72hpi_adult_data.filt_1 <- subset(Short_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample1_Sample1")
Short_read_VIC01_72hpi_adult_data.filt_2 <- subset(Short_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample2_Sample2")
Short_read_VIC01_72hpi_adult_data.filt_3 <- subset(Short_read_VIC01_72hpi_adult_data.filt, subset = Sample =="Sample3_Sample3")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_VIC01_72hpi_adult_data.filt_1 <- NormalizeData(Short_read_VIC01_72hpi_adult_data.filt_1)
Short_read_VIC01_72hpi_adult_data.filt_2 <- NormalizeData(Short_read_VIC01_72hpi_adult_data.filt_2)
Short_read_VIC01_72hpi_adult_data.filt_3 <- NormalizeData(Short_read_VIC01_72hpi_adult_data.filt_3)



Short_read_VIC01_72hpi_adult_data.filt_1 <- FindVariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_1, selection.method = "vst")
Short_read_VIC01_72hpi_adult_data.filt_2 <- FindVariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_2, selection.method = "vst")
Short_read_VIC01_72hpi_adult_data.filt_3 <- FindVariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_3, selection.method = "vst")


Short_read_VIC01_72hpi_adult_data.filt_1 <- ScaleData(Short_read_VIC01_72hpi_adult_data.filt_1, features = rownames(Short_read_VIC01_72hpi_adult_data.filt_1))
Short_read_VIC01_72hpi_adult_data.filt_2 <- ScaleData(Short_read_VIC01_72hpi_adult_data.filt_2, features = rownames(Short_read_VIC01_72hpi_adult_data.filt_2))
Short_read_VIC01_72hpi_adult_data.filt_3 <- ScaleData(Short_read_VIC01_72hpi_adult_data.filt_3, features = rownames(Short_read_VIC01_72hpi_adult_data.filt_3))



Short_read_VIC01_72hpi_adult_data.filt_1 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_1, features = VariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_72hpi_adult_data.filt_2 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_2, features = VariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_72hpi_adult_data.filt_3 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_3, features = VariableFeatures(Short_read_VIC01_72hpi_adult_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)


Short_read_VIC01_72hpi_adult_data.filt_1 <- CellCycleScoring(Short_read_VIC01_72hpi_adult_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_72hpi_adult_data.filt_2 <- CellCycleScoring(Short_read_VIC01_72hpi_adult_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_72hpi_adult_data.filt_3 <- CellCycleScoring(Short_read_VIC01_72hpi_adult_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Short_read_VIC01_72hpi_adult_data.filt_1$CC.Difference <- Short_read_VIC01_72hpi_adult_data.filt_1$S.Score - Short_read_VIC01_72hpi_adult_data.filt_1$G2M.Score
Short_read_VIC01_72hpi_adult_data.filt_2$CC.Difference <- Short_read_VIC01_72hpi_adult_data.filt_2$S.Score - Short_read_VIC01_72hpi_adult_data.filt_2$G2M.Score
Short_read_VIC01_72hpi_adult_data.filt_3$CC.Difference <- Short_read_VIC01_72hpi_adult_data.filt_3$S.Score - Short_read_VIC01_72hpi_adult_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_VIC01_72hpi_adult_data.filt_1 <- SCTransform(Short_read_VIC01_72hpi_adult_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_72hpi_adult_data.filt_2 <- SCTransform(Short_read_VIC01_72hpi_adult_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_72hpi_adult_data.filt_3 <- SCTransform(Short_read_VIC01_72hpi_adult_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_VIC01_72hpi_adult_data.filt_1 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_72hpi_adult_data.filt_2 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_72hpi_adult_data.filt_3 <- RunPCA(Short_read_VIC01_72hpi_adult_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Short_read_VIC01_72hpi_adult_data.filt_1)
DimPlot(Short_read_VIC01_72hpi_adult_data.filt_2)
DimPlot(Short_read_VIC01_72hpi_adult_data.filt_3)


# store mitochondrial percentage in object meta data
Long_read_UK_72hpi_child <- PercentageFeatureSet(Long_read_UK_72hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Long_read_UK_72hpi_child_f <- PercentageFeatureSet(Long_read_UK_72hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_UK_72hpi_child <- PercentageFeatureSet(Long_read_UK_72hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_UK_72hpi_child <- PercentageFeatureSet(Long_read_UK_72hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_UK_72hpi_child <- SCTransform(Long_read_UK_72hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_UK_72hpi_child_f <- SCTransform(Long_read_UK_72hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_UK_72hpi_child_data <- subset(Long_read_UK_72hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_UK_72hpi_child_selected_c <- WhichCells(Long_read_UK_72hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_UK_72hpi_child_selected_f <- rownames(Long_read_UK_72hpi_child_data)[Matrix::rowSums(Long_read_UK_72hpi_child_data) > 3]

Long_read_UK_72hpi_child_data.filt <- subset(Long_read_UK_72hpi_child_data, features = Long_read_UK_72hpi_child_selected_f, cells = Long_read_UK_72hpi_child_selected_c)
dim(Long_read_UK_72hpi_child_data.filt)




#mito/ribo filtering
Long_read_UK_72hpi_child_selected_mito <- WhichCells(Long_read_UK_72hpi_child_data.filt, expression = percent_mito < 20)
Long_read_UK_72hpi_child_selected_ribo <- WhichCells(Long_read_UK_72hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_UK_72hpi_child_data.filt <- subset(Long_read_UK_72hpi_child_data.filt, cells = Long_read_UK_72hpi_child_selected_mito)
Long_read_UK_72hpi_child_data.filt <- subset(Long_read_UK_72hpi_child_data.filt, cells = Long_read_UK_72hpi_child_selected_ribo)

dim(Long_read_UK_72hpi_child_data.filt)

table(Long_read_UK_72hpi_child_data.filt$orig.ident)


Long_read_UK_72hpi_child_data.filt_1 <- subset(Long_read_UK_72hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Long_read_UK_72hpi_child_data.filt_2 <- subset(Long_read_UK_72hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Long_read_UK_72hpi_child_data.filt_3 <- subset(Long_read_UK_72hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_UK_72hpi_child_data.filt_1 <- NormalizeData(Long_read_UK_72hpi_child_data.filt_1)
Long_read_UK_72hpi_child_data.filt_2 <- NormalizeData(Long_read_UK_72hpi_child_data.filt_2)
Long_read_UK_72hpi_child_data.filt_3 <- NormalizeData(Long_read_UK_72hpi_child_data.filt_3)



Long_read_UK_72hpi_child_data.filt_1 <- FindVariableFeatures(Long_read_UK_72hpi_child_data.filt_1, selection.method = "vst")
Long_read_UK_72hpi_child_data.filt_2 <- FindVariableFeatures(Long_read_UK_72hpi_child_data.filt_2, selection.method = "vst")
Long_read_UK_72hpi_child_data.filt_3 <- FindVariableFeatures(Long_read_UK_72hpi_child_data.filt_3, selection.method = "vst")


Long_read_UK_72hpi_child_data.filt_1 <- ScaleData(Long_read_UK_72hpi_child_data.filt_1, features = rownames(Long_read_UK_72hpi_child_data.filt_1))
Long_read_UK_72hpi_child_data.filt_2 <- ScaleData(Long_read_UK_72hpi_child_data.filt_2, features = rownames(Long_read_UK_72hpi_child_data.filt_2))
Long_read_UK_72hpi_child_data.filt_3 <- ScaleData(Long_read_UK_72hpi_child_data.filt_3, features = rownames(Long_read_UK_72hpi_child_data.filt_3))



Long_read_UK_72hpi_child_data.filt_1 <- RunPCA(Long_read_UK_72hpi_child_data.filt_1, features = VariableFeatures(Long_read_UK_72hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_UK_72hpi_child_data.filt_2 <- RunPCA(Long_read_UK_72hpi_child_data.filt_2, features = VariableFeatures(Long_read_UK_72hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_UK_72hpi_child_data.filt_3 <- RunPCA(Long_read_UK_72hpi_child_data.filt_3, features = VariableFeatures(Long_read_UK_72hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)





Long_read_UK_72hpi_child_data.filt_1 <- CellCycleScoring(Long_read_UK_72hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_UK_72hpi_child_data.filt_2 <- CellCycleScoring(Long_read_UK_72hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_UK_72hpi_child_data.filt_3 <- CellCycleScoring(Long_read_UK_72hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



#Alternate Workflow regression
Long_read_UK_72hpi_child_data.filt_1$CC.Difference <- Long_read_UK_72hpi_child_data.filt_1$S.Score - Long_read_UK_72hpi_child_data.filt_1$G2M.Score
Long_read_UK_72hpi_child_data.filt_2$CC.Difference <- Long_read_UK_72hpi_child_data.filt_2$S.Score - Long_read_UK_72hpi_child_data.filt_2$G2M.Score
Long_read_UK_72hpi_child_data.filt_3$CC.Difference <- Long_read_UK_72hpi_child_data.filt_3$S.Score - Long_read_UK_72hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_UK_72hpi_child_data.filt_1 <- SCTransform(Long_read_UK_72hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_UK_72hpi_child_data.filt_2 <- SCTransform(Long_read_UK_72hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_UK_72hpi_child_data.filt_3 <- SCTransform(Long_read_UK_72hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_UK_72hpi_child_data.filt_1 <- RunPCA(Long_read_UK_72hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_UK_72hpi_child_data.filt_2 <- RunPCA(Long_read_UK_72hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_UK_72hpi_child_data.filt_3 <- RunPCA(Long_read_UK_72hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Long_read_UK_72hpi_child_data.filt_1)
DimPlot(Long_read_UK_72hpi_child_data.filt_2)
DimPlot(Long_read_UK_72hpi_child_data.filt_3)


# store mitochondrial percentage in object meta data
Long_read_uninfected_child <- PercentageFeatureSet(Long_read_uninfected_child, pattern = "^MT-", col.name = "percent_mito")
#Long_read_uninfected_child_f <- PercentageFeatureSet(Long_read_uninfected_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_uninfected_child <- PercentageFeatureSet(Long_read_uninfected_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_uninfected_child <- PercentageFeatureSet(Long_read_uninfected_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_uninfected_child <- SCTransform(Long_read_uninfected_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_uninfected_child_f <- SCTransform(Long_read_uninfected_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_uninfected_child_data <- subset(Long_read_uninfected_child, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_uninfected_child_selected_c <- WhichCells(Long_read_uninfected_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_uninfected_child_selected_f <- rownames(Long_read_uninfected_child_data)[Matrix::rowSums(Long_read_uninfected_child_data) > 3]

Long_read_uninfected_child_data.filt <- subset(Long_read_uninfected_child_data, features = Long_read_uninfected_child_selected_f, cells = Long_read_uninfected_child_selected_c)
dim(Long_read_uninfected_child_data.filt)




#mito/ribo filtering
Long_read_uninfected_child_selected_mito <- WhichCells(Long_read_uninfected_child_data.filt, expression = percent_mito < 20)
Long_read_uninfected_child_selected_ribo <- WhichCells(Long_read_uninfected_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_uninfected_child_data.filt <- subset(Long_read_uninfected_child_data.filt, cells = Long_read_uninfected_child_selected_mito)
Long_read_uninfected_child_data.filt <- subset(Long_read_uninfected_child_data.filt, cells = Long_read_uninfected_child_selected_ribo)

dim(Long_read_uninfected_child_data.filt)

table(Long_read_uninfected_child_data.filt$orig.ident)


Long_read_uninfected_child_data.filt_1 <- subset(Long_read_uninfected_child_data.filt, subset = Sample =="Sample4_Sample4")
Long_read_uninfected_child_data.filt_2 <- subset(Long_read_uninfected_child_data.filt, subset = Sample =="Sample5_Sample5")
Long_read_uninfected_child_data.filt_3 <- subset(Long_read_uninfected_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_uninfected_child_data.filt_1 <- NormalizeData(Long_read_uninfected_child_data.filt_1)
Long_read_uninfected_child_data.filt_2 <- NormalizeData(Long_read_uninfected_child_data.filt_2)
Long_read_uninfected_child_data.filt_3 <- NormalizeData(Long_read_uninfected_child_data.filt_3)



Long_read_uninfected_child_data.filt_1 <- FindVariableFeatures(Long_read_uninfected_child_data.filt_1, selection.method = "vst")
Long_read_uninfected_child_data.filt_2 <- FindVariableFeatures(Long_read_uninfected_child_data.filt_2, selection.method = "vst")
Long_read_uninfected_child_data.filt_3 <- FindVariableFeatures(Long_read_uninfected_child_data.filt_3, selection.method = "vst")


Long_read_uninfected_child_data.filt_1 <- ScaleData(Long_read_uninfected_child_data.filt_1, features = rownames(Long_read_uninfected_child_data.filt_1))
Long_read_uninfected_child_data.filt_2 <- ScaleData(Long_read_uninfected_child_data.filt_2, features = rownames(Long_read_uninfected_child_data.filt_2))
Long_read_uninfected_child_data.filt_3 <- ScaleData(Long_read_uninfected_child_data.filt_3, features = rownames(Long_read_uninfected_child_data.filt_3))



Long_read_uninfected_child_data.filt_1 <- RunPCA(Long_read_uninfected_child_data.filt_1, features = VariableFeatures(Long_read_uninfected_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_uninfected_child_data.filt_2 <- RunPCA(Long_read_uninfected_child_data.filt_2, features = VariableFeatures(Long_read_uninfected_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_uninfected_child_data.filt_3 <- RunPCA(Long_read_uninfected_child_data.filt_3, features = VariableFeatures(Long_read_uninfected_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)




Long_read_uninfected_child_data.filt_1 <- CellCycleScoring(Long_read_uninfected_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_uninfected_child_data.filt_2 <- CellCycleScoring(Long_read_uninfected_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_uninfected_child_data.filt_3 <- CellCycleScoring(Long_read_uninfected_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Long_read_uninfected_child_data.filt_1$CC.Difference <- Long_read_uninfected_child_data.filt_1$S.Score - Long_read_uninfected_child_data.filt_1$G2M.Score
Long_read_uninfected_child_data.filt_2$CC.Difference <- Long_read_uninfected_child_data.filt_2$S.Score - Long_read_uninfected_child_data.filt_2$G2M.Score
Long_read_uninfected_child_data.filt_3$CC.Difference <- Long_read_uninfected_child_data.filt_3$S.Score - Long_read_uninfected_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_uninfected_child_data.filt_1 <- SCTransform(Long_read_uninfected_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_uninfected_child_data.filt_2 <- SCTransform(Long_read_uninfected_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_uninfected_child_data.filt_3 <- SCTransform(Long_read_uninfected_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_uninfected_child_data.filt_1 <- RunPCA(Long_read_uninfected_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_uninfected_child_data.filt_2 <- RunPCA(Long_read_uninfected_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_uninfected_child_data.filt_3 <- RunPCA(Long_read_uninfected_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Long_read_uninfected_child_data.filt_1)
DimPlot(Long_read_uninfected_child_data.filt_2)
DimPlot(Long_read_uninfected_child_data.filt_3)


# store mitochondrial percentage in object meta data
Long_read_VIC01_48hpi_child <- PercentageFeatureSet(Long_read_VIC01_48hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Long_read_VIC01_48hpi_child_f <- PercentageFeatureSet(Long_read_VIC01_48hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_VIC01_48hpi_child <- PercentageFeatureSet(Long_read_VIC01_48hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_VIC01_48hpi_child <- PercentageFeatureSet(Long_read_VIC01_48hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_VIC01_48hpi_child <- SCTransform(Long_read_VIC01_48hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_VIC01_48hpi_child_f <- SCTransform(Long_read_VIC01_48hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_VIC01_48hpi_child_data <- subset(Long_read_VIC01_48hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_VIC01_48hpi_child_selected_c <- WhichCells(Long_read_VIC01_48hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_VIC01_48hpi_child_selected_f <- rownames(Long_read_VIC01_48hpi_child_data)[Matrix::rowSums(Long_read_VIC01_48hpi_child_data) > 3]

Long_read_VIC01_48hpi_child_data.filt <- subset(Long_read_VIC01_48hpi_child_data, features = Long_read_VIC01_48hpi_child_selected_f, cells = Long_read_VIC01_48hpi_child_selected_c)
dim(Long_read_VIC01_48hpi_child_data.filt)




#mito/ribo filtering
Long_read_VIC01_48hpi_child_selected_mito <- WhichCells(Long_read_VIC01_48hpi_child_data.filt, expression = percent_mito < 20)
Long_read_VIC01_48hpi_child_selected_ribo <- WhichCells(Long_read_VIC01_48hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_VIC01_48hpi_child_data.filt <- subset(Long_read_VIC01_48hpi_child_data.filt, cells = Long_read_VIC01_48hpi_child_selected_mito)
Long_read_VIC01_48hpi_child_data.filt <- subset(Long_read_VIC01_48hpi_child_data.filt, cells = Long_read_VIC01_48hpi_child_selected_ribo)

dim(Long_read_VIC01_48hpi_child_data.filt)

table(Long_read_VIC01_48hpi_child_data.filt$orig.ident)


Long_read_VIC01_48hpi_child_data.filt_1 <- subset(Long_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Long_read_VIC01_48hpi_child_data.filt_2 <- subset(Long_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Long_read_VIC01_48hpi_child_data.filt_3 <- subset(Long_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_VIC01_48hpi_child_data.filt_1 <- NormalizeData(Long_read_VIC01_48hpi_child_data.filt_1)
Long_read_VIC01_48hpi_child_data.filt_2 <- NormalizeData(Long_read_VIC01_48hpi_child_data.filt_2)
Long_read_VIC01_48hpi_child_data.filt_3 <- NormalizeData(Long_read_VIC01_48hpi_child_data.filt_3)



Long_read_VIC01_48hpi_child_data.filt_1 <- FindVariableFeatures(Long_read_VIC01_48hpi_child_data.filt_1, selection.method = "vst")
Long_read_VIC01_48hpi_child_data.filt_2 <- FindVariableFeatures(Long_read_VIC01_48hpi_child_data.filt_2, selection.method = "vst")
Long_read_VIC01_48hpi_child_data.filt_3 <- FindVariableFeatures(Long_read_VIC01_48hpi_child_data.filt_3, selection.method = "vst")


Long_read_VIC01_48hpi_child_data.filt_1 <- ScaleData(Long_read_VIC01_48hpi_child_data.filt_1, features = rownames(Long_read_VIC01_48hpi_child_data.filt_1))
Long_read_VIC01_48hpi_child_data.filt_2 <- ScaleData(Long_read_VIC01_48hpi_child_data.filt_2, features = rownames(Long_read_VIC01_48hpi_child_data.filt_2))
Long_read_VIC01_48hpi_child_data.filt_3 <- ScaleData(Long_read_VIC01_48hpi_child_data.filt_3, features = rownames(Long_read_VIC01_48hpi_child_data.filt_3))



Long_read_VIC01_48hpi_child_data.filt_1 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_1, features = VariableFeatures(Long_read_VIC01_48hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_48hpi_child_data.filt_2 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_2, features = VariableFeatures(Long_read_VIC01_48hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_48hpi_child_data.filt_3 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_3, features = VariableFeatures(Long_read_VIC01_48hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)




Long_read_VIC01_48hpi_child_data.filt_1 <- CellCycleScoring(Long_read_VIC01_48hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_48hpi_child_data.filt_2 <- CellCycleScoring(Long_read_VIC01_48hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_48hpi_child_data.filt_3 <- CellCycleScoring(Long_read_VIC01_48hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Long_read_VIC01_48hpi_child_data.filt_1$CC.Difference <- Long_read_VIC01_48hpi_child_data.filt_1$S.Score - Long_read_VIC01_48hpi_child_data.filt_1$G2M.Score
Long_read_VIC01_48hpi_child_data.filt_2$CC.Difference <- Long_read_VIC01_48hpi_child_data.filt_2$S.Score - Long_read_VIC01_48hpi_child_data.filt_2$G2M.Score
Long_read_VIC01_48hpi_child_data.filt_3$CC.Difference <- Long_read_VIC01_48hpi_child_data.filt_3$S.Score - Long_read_VIC01_48hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_VIC01_48hpi_child_data.filt_1 <- SCTransform(Long_read_VIC01_48hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_48hpi_child_data.filt_2 <- SCTransform(Long_read_VIC01_48hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_48hpi_child_data.filt_3 <- SCTransform(Long_read_VIC01_48hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_VIC01_48hpi_child_data.filt_1 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_48hpi_child_data.filt_2 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_48hpi_child_data.filt_3 <- RunPCA(Long_read_VIC01_48hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Long_read_VIC01_48hpi_child_data.filt_1)
DimPlot(Long_read_VIC01_48hpi_child_data.filt_2)
DimPlot(Long_read_VIC01_48hpi_child_data.filt_3)


# store mitochondrial percentage in object meta data
Long_read_VIC01_72hpi_child <- PercentageFeatureSet(Long_read_VIC01_72hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Long_read_VIC01_72hpi_child_f <- PercentageFeatureSet(Long_read_VIC01_72hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Long_read_VIC01_72hpi_child <- PercentageFeatureSet(Long_read_VIC01_72hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Long_read_VIC01_72hpi_child <- PercentageFeatureSet(Long_read_VIC01_72hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Long_read_VIC01_72hpi_child <- SCTransform(Long_read_VIC01_72hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Long_read_VIC01_72hpi_child_f <- SCTransform(Long_read_VIC01_72hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Long_read_VIC01_72hpi_child_data <- subset(Long_read_VIC01_72hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Long_read_VIC01_72hpi_child_selected_c <- WhichCells(Long_read_VIC01_72hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Long_read_VIC01_72hpi_child_selected_f <- rownames(Long_read_VIC01_72hpi_child_data)[Matrix::rowSums(Long_read_VIC01_72hpi_child_data) > 3]

Long_read_VIC01_72hpi_child_data.filt <- subset(Long_read_VIC01_72hpi_child_data, features = Long_read_VIC01_72hpi_child_selected_f, cells = Long_read_VIC01_72hpi_child_selected_c)
dim(Long_read_VIC01_72hpi_child_data.filt)




#mito/ribo filtering
Long_read_VIC01_72hpi_child_selected_mito <- WhichCells(Long_read_VIC01_72hpi_child_data.filt, expression = percent_mito < 20)
Long_read_VIC01_72hpi_child_selected_ribo <- WhichCells(Long_read_VIC01_72hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Long_read_VIC01_72hpi_child_data.filt <- subset(Long_read_VIC01_72hpi_child_data.filt, cells = Long_read_VIC01_72hpi_child_selected_mito)
Long_read_VIC01_72hpi_child_data.filt <- subset(Long_read_VIC01_72hpi_child_data.filt, cells = Long_read_VIC01_72hpi_child_selected_ribo)

dim(Long_read_VIC01_72hpi_child_data.filt)

table(Long_read_VIC01_72hpi_child_data.filt$orig.ident)


Long_read_VIC01_72hpi_child_data.filt_1 <- subset(Long_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Long_read_VIC01_72hpi_child_data.filt_2 <- subset(Long_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Long_read_VIC01_72hpi_child_data.filt_3 <- subset(Long_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Long_read_VIC01_72hpi_child_data.filt_1 <- NormalizeData(Long_read_VIC01_72hpi_child_data.filt_1)
Long_read_VIC01_72hpi_child_data.filt_2 <- NormalizeData(Long_read_VIC01_72hpi_child_data.filt_2)
Long_read_VIC01_72hpi_child_data.filt_3 <- NormalizeData(Long_read_VIC01_72hpi_child_data.filt_3)



Long_read_VIC01_72hpi_child_data.filt_1 <- FindVariableFeatures(Long_read_VIC01_72hpi_child_data.filt_1, selection.method = "vst")
Long_read_VIC01_72hpi_child_data.filt_2 <- FindVariableFeatures(Long_read_VIC01_72hpi_child_data.filt_2, selection.method = "vst")
Long_read_VIC01_72hpi_child_data.filt_3 <- FindVariableFeatures(Long_read_VIC01_72hpi_child_data.filt_3, selection.method = "vst")


Long_read_VIC01_72hpi_child_data.filt_1 <- ScaleData(Long_read_VIC01_72hpi_child_data.filt_1, features = rownames(Long_read_VIC01_72hpi_child_data.filt_1))
Long_read_VIC01_72hpi_child_data.filt_2 <- ScaleData(Long_read_VIC01_72hpi_child_data.filt_2, features = rownames(Long_read_VIC01_72hpi_child_data.filt_2))
Long_read_VIC01_72hpi_child_data.filt_3 <- ScaleData(Long_read_VIC01_72hpi_child_data.filt_3, features = rownames(Long_read_VIC01_72hpi_child_data.filt_3))



Long_read_VIC01_72hpi_child_data.filt_1 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_1, features = VariableFeatures(Long_read_VIC01_72hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_72hpi_child_data.filt_2 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_2, features = VariableFeatures(Long_read_VIC01_72hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Long_read_VIC01_72hpi_child_data.filt_3 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_3, features = VariableFeatures(Long_read_VIC01_72hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)





Long_read_VIC01_72hpi_child_data.filt_1 <- CellCycleScoring(Long_read_VIC01_72hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_72hpi_child_data.filt_2 <- CellCycleScoring(Long_read_VIC01_72hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Long_read_VIC01_72hpi_child_data.filt_3 <- CellCycleScoring(Long_read_VIC01_72hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Long_read_VIC01_72hpi_child_data.filt_1$CC.Difference <- Long_read_VIC01_72hpi_child_data.filt_1$S.Score - Long_read_VIC01_72hpi_child_data.filt_1$G2M.Score
Long_read_VIC01_72hpi_child_data.filt_2$CC.Difference <- Long_read_VIC01_72hpi_child_data.filt_2$S.Score - Long_read_VIC01_72hpi_child_data.filt_2$G2M.Score
Long_read_VIC01_72hpi_child_data.filt_3$CC.Difference <- Long_read_VIC01_72hpi_child_data.filt_3$S.Score - Long_read_VIC01_72hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Long_read_VIC01_72hpi_child_data.filt_1 <- SCTransform(Long_read_VIC01_72hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_72hpi_child_data.filt_2 <- SCTransform(Long_read_VIC01_72hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Long_read_VIC01_72hpi_child_data.filt_3 <- SCTransform(Long_read_VIC01_72hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Long_read_VIC01_72hpi_child_data.filt_1 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_72hpi_child_data.filt_2 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Long_read_VIC01_72hpi_child_data.filt_3 <- RunPCA(Long_read_VIC01_72hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Long_read_VIC01_72hpi_child_data.filt_1)
DimPlot(Long_read_VIC01_72hpi_child_data.filt_2)
DimPlot(Long_read_VIC01_72hpi_child_data.filt_3)



# store mitochondrial percentage in object meta data
Short_read_UK_72hpi_child <- PercentageFeatureSet(Short_read_UK_72hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Short_read_UK_72hpi_child_f <- PercentageFeatureSet(Short_read_UK_72hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_UK_72hpi_child <- PercentageFeatureSet(Short_read_UK_72hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_UK_72hpi_child <- PercentageFeatureSet(Short_read_UK_72hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_UK_72hpi_child <- SCTransform(Short_read_UK_72hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_UK_72hpi_child_f <- SCTransform(Short_read_UK_72hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_UK_72hpi_child_data <- subset(Short_read_UK_72hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_UK_72hpi_child_selected_c <- WhichCells(Short_read_UK_72hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_UK_72hpi_child_selected_f <- rownames(Short_read_UK_72hpi_child_data)[Matrix::rowSums(Short_read_UK_72hpi_child_data) > 3]

Short_read_UK_72hpi_child_data.filt <- subset(Short_read_UK_72hpi_child_data, features = Short_read_UK_72hpi_child_selected_f, cells = Short_read_UK_72hpi_child_selected_c)
dim(Short_read_UK_72hpi_child_data.filt)




#mito/ribo filtering
Short_read_UK_72hpi_child_selected_mito <- WhichCells(Short_read_UK_72hpi_child_data.filt, expression = percent_mito < 20)
Short_read_UK_72hpi_child_selected_ribo <- WhichCells(Short_read_UK_72hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_UK_72hpi_child_data.filt <- subset(Short_read_UK_72hpi_child_data.filt, cells = Short_read_UK_72hpi_child_selected_mito)
Short_read_UK_72hpi_child_data.filt <- subset(Short_read_UK_72hpi_child_data.filt, cells = Short_read_UK_72hpi_child_selected_ribo)

dim(Short_read_UK_72hpi_child_data.filt)

table(Short_read_UK_72hpi_child_data.filt$orig.ident)


Short_read_UK_72hpi_child_data.filt_1 <- subset(Short_read_UK_72hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Short_read_UK_72hpi_child_data.filt_2 <- subset(Short_read_UK_72hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Short_read_UK_72hpi_child_data.filt_3 <- subset(Short_read_UK_72hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_UK_72hpi_child_data.filt_1 <- NormalizeData(Short_read_UK_72hpi_child_data.filt_1)
Short_read_UK_72hpi_child_data.filt_2 <- NormalizeData(Short_read_UK_72hpi_child_data.filt_2)
Short_read_UK_72hpi_child_data.filt_3 <- NormalizeData(Short_read_UK_72hpi_child_data.filt_3)



Short_read_UK_72hpi_child_data.filt_1 <- FindVariableFeatures(Short_read_UK_72hpi_child_data.filt_1, selection.method = "vst")
Short_read_UK_72hpi_child_data.filt_2 <- FindVariableFeatures(Short_read_UK_72hpi_child_data.filt_2, selection.method = "vst")
Short_read_UK_72hpi_child_data.filt_3 <- FindVariableFeatures(Short_read_UK_72hpi_child_data.filt_3, selection.method = "vst")


Short_read_UK_72hpi_child_data.filt_1 <- ScaleData(Short_read_UK_72hpi_child_data.filt_1, features = rownames(Short_read_UK_72hpi_child_data.filt_1))
Short_read_UK_72hpi_child_data.filt_2 <- ScaleData(Short_read_UK_72hpi_child_data.filt_2, features = rownames(Short_read_UK_72hpi_child_data.filt_2))
Short_read_UK_72hpi_child_data.filt_3 <- ScaleData(Short_read_UK_72hpi_child_data.filt_3, features = rownames(Short_read_UK_72hpi_child_data.filt_3))



Short_read_UK_72hpi_child_data.filt_1 <- RunPCA(Short_read_UK_72hpi_child_data.filt_1, features = VariableFeatures(Short_read_UK_72hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_UK_72hpi_child_data.filt_2 <- RunPCA(Short_read_UK_72hpi_child_data.filt_2, features = VariableFeatures(Short_read_UK_72hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_UK_72hpi_child_data.filt_3 <- RunPCA(Short_read_UK_72hpi_child_data.filt_3, features = VariableFeatures(Short_read_UK_72hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)




Short_read_UK_72hpi_child_data.filt_1 <- CellCycleScoring(Short_read_UK_72hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_UK_72hpi_child_data.filt_2 <- CellCycleScoring(Short_read_UK_72hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_UK_72hpi_child_data.filt_3 <- CellCycleScoring(Short_read_UK_72hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#Alternate Workflow regression
Short_read_UK_72hpi_child_data.filt_1$CC.Difference <- Short_read_UK_72hpi_child_data.filt_1$S.Score - Short_read_UK_72hpi_child_data.filt_1$G2M.Score
Short_read_UK_72hpi_child_data.filt_2$CC.Difference <- Short_read_UK_72hpi_child_data.filt_2$S.Score - Short_read_UK_72hpi_child_data.filt_2$G2M.Score
Short_read_UK_72hpi_child_data.filt_3$CC.Difference <- Short_read_UK_72hpi_child_data.filt_3$S.Score - Short_read_UK_72hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_UK_72hpi_child_data.filt_1 <- SCTransform(Short_read_UK_72hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_UK_72hpi_child_data.filt_2 <- SCTransform(Short_read_UK_72hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_UK_72hpi_child_data.filt_3 <- SCTransform(Short_read_UK_72hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_UK_72hpi_child_data.filt_1 <- RunPCA(Short_read_UK_72hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_UK_72hpi_child_data.filt_2 <- RunPCA(Short_read_UK_72hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_UK_72hpi_child_data.filt_3 <- RunPCA(Short_read_UK_72hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_UK_72hpi_child_data.filt_1)
DimPlot(Short_read_UK_72hpi_child_data.filt_2)
DimPlot(Short_read_UK_72hpi_child_data.filt_3)


# store mitochondrial percentage in object meta data
Short_read_uninfected_child<- PercentageFeatureSet(Short_read_uninfected_child, pattern = "^MT-", col.name = "percent_mito")
#Short_read_uninfected_child_f <- PercentageFeatureSet(Short_read_uninfected_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_uninfected_child<- PercentageFeatureSet(Short_read_uninfected_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_uninfected_child<- PercentageFeatureSet(Short_read_uninfected_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_uninfected_child<- SCTransform(Short_read_uninfected_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_uninfected_child_f <- SCTransform(Short_read_uninfected_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_uninfected_child_data <- subset(Short_read_uninfected_child, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_uninfected_child_selected_c <- WhichCells(Short_read_uninfected_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_uninfected_child_selected_f <- rownames(Short_read_uninfected_child_data)[Matrix::rowSums(Short_read_uninfected_child_data) > 3]

Short_read_uninfected_child_data.filt <- subset(Short_read_uninfected_child_data, features = Short_read_uninfected_child_selected_f, cells = Short_read_uninfected_child_selected_c)
dim(Short_read_uninfected_child_data.filt)




#mito/ribo filtering
Short_read_uninfected_child_selected_mito <- WhichCells(Short_read_uninfected_child_data.filt, expression = percent_mito < 20)
Short_read_uninfected_child_selected_ribo <- WhichCells(Short_read_uninfected_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_uninfected_child_data.filt <- subset(Short_read_uninfected_child_data.filt, cells = Short_read_uninfected_child_selected_mito)
Short_read_uninfected_child_data.filt <- subset(Short_read_uninfected_child_data.filt, cells = Short_read_uninfected_child_selected_ribo)

dim(Short_read_uninfected_child_data.filt)

table(Short_read_uninfected_child_data.filt$orig.ident)


Short_read_uninfected_child_data.filt_1 <- subset(Short_read_uninfected_child_data.filt, subset = Sample =="Sample4_Sample4")
Short_read_uninfected_child_data.filt_2 <- subset(Short_read_uninfected_child_data.filt, subset = Sample =="Sample5_Sample5")
Short_read_uninfected_child_data.filt_3 <- subset(Short_read_uninfected_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_uninfected_child_data.filt_1 <- NormalizeData(Short_read_uninfected_child_data.filt_1)
Short_read_uninfected_child_data.filt_2 <- NormalizeData(Short_read_uninfected_child_data.filt_2)
Short_read_uninfected_child_data.filt_3 <- NormalizeData(Short_read_uninfected_child_data.filt_3)



Short_read_uninfected_child_data.filt_1 <- FindVariableFeatures(Short_read_uninfected_child_data.filt_1, selection.method = "vst")
Short_read_uninfected_child_data.filt_2 <- FindVariableFeatures(Short_read_uninfected_child_data.filt_2, selection.method = "vst")
Short_read_uninfected_child_data.filt_3 <- FindVariableFeatures(Short_read_uninfected_child_data.filt_3, selection.method = "vst")


Short_read_uninfected_child_data.filt_1 <- ScaleData(Short_read_uninfected_child_data.filt_1, features = rownames(Short_read_uninfected_child_data.filt_1))
Short_read_uninfected_child_data.filt_2 <- ScaleData(Short_read_uninfected_child_data.filt_2, features = rownames(Short_read_uninfected_child_data.filt_2))
Short_read_uninfected_child_data.filt_3 <- ScaleData(Short_read_uninfected_child_data.filt_3, features = rownames(Short_read_uninfected_child_data.filt_3))



Short_read_uninfected_child_data.filt_1 <- RunPCA(Short_read_uninfected_child_data.filt_1, features = VariableFeatures(Short_read_uninfected_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_uninfected_child_data.filt_2 <- RunPCA(Short_read_uninfected_child_data.filt_2, features = VariableFeatures(Short_read_uninfected_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_uninfected_child_data.filt_3 <- RunPCA(Short_read_uninfected_child_data.filt_3, features = VariableFeatures(Short_read_uninfected_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)



Short_read_uninfected_child_data.filt_1 <- CellCycleScoring(Short_read_uninfected_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_uninfected_child_data.filt_2 <- CellCycleScoring(Short_read_uninfected_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_uninfected_child_data.filt_3 <- CellCycleScoring(Short_read_uninfected_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Short_read_uninfected_child_data.filt_1$CC.Difference <- Short_read_uninfected_child_data.filt_1$S.Score - Short_read_uninfected_child_data.filt_1$G2M.Score
Short_read_uninfected_child_data.filt_2$CC.Difference <- Short_read_uninfected_child_data.filt_2$S.Score - Short_read_uninfected_child_data.filt_2$G2M.Score
Short_read_uninfected_child_data.filt_3$CC.Difference <- Short_read_uninfected_child_data.filt_3$S.Score - Short_read_uninfected_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_uninfected_child_data.filt_1 <- SCTransform(Short_read_uninfected_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_uninfected_child_data.filt_2 <- SCTransform(Short_read_uninfected_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_uninfected_child_data.filt_3 <- SCTransform(Short_read_uninfected_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_uninfected_child_data.filt_1 <- RunPCA(Short_read_uninfected_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_uninfected_child_data.filt_2 <- RunPCA(Short_read_uninfected_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_uninfected_child_data.filt_3 <- RunPCA(Short_read_uninfected_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_uninfected_child_data.filt_1)
DimPlot(Short_read_uninfected_child_data.filt_2)
DimPlot(Short_read_uninfected_child_data.filt_3)


# store mitochondrial percentage in object meta data
Short_read_VIC01_48hpi_child<- PercentageFeatureSet(Short_read_VIC01_48hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Short_read_VIC01_48hpi_child_f <- PercentageFeatureSet(Short_read_VIC01_48hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_VIC01_48hpi_child<- PercentageFeatureSet(Short_read_VIC01_48hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_VIC01_48hpi_child<- PercentageFeatureSet(Short_read_VIC01_48hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_VIC01_48hpi_child<- SCTransform(Short_read_VIC01_48hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_VIC01_48hpi_child_f <- SCTransform(Short_read_VIC01_48hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_VIC01_48hpi_child_data <- subset(Short_read_VIC01_48hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_VIC01_48hpi_child_selected_c <- WhichCells(Short_read_VIC01_48hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_VIC01_48hpi_child_selected_f <- rownames(Short_read_VIC01_48hpi_child_data)[Matrix::rowSums(Short_read_VIC01_48hpi_child_data) > 3]

Short_read_VIC01_48hpi_child_data.filt <- subset(Short_read_VIC01_48hpi_child_data, features = Short_read_VIC01_48hpi_child_selected_f, cells = Short_read_VIC01_48hpi_child_selected_c)
dim(Short_read_VIC01_48hpi_child_data.filt)




#mito/ribo filtering
Short_read_VIC01_48hpi_child_selected_mito <- WhichCells(Short_read_VIC01_48hpi_child_data.filt, expression = percent_mito < 20)
Short_read_VIC01_48hpi_child_selected_ribo <- WhichCells(Short_read_VIC01_48hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_VIC01_48hpi_child_data.filt <- subset(Short_read_VIC01_48hpi_child_data.filt, cells = Short_read_VIC01_48hpi_child_selected_mito)
Short_read_VIC01_48hpi_child_data.filt <- subset(Short_read_VIC01_48hpi_child_data.filt, cells = Short_read_VIC01_48hpi_child_selected_ribo)

dim(Short_read_VIC01_48hpi_child_data.filt)

table(Short_read_VIC01_48hpi_child_data.filt$orig.ident)


Short_read_VIC01_48hpi_child_data.filt_1 <- subset(Short_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Short_read_VIC01_48hpi_child_data.filt_2 <- subset(Short_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Short_read_VIC01_48hpi_child_data.filt_3 <- subset(Short_read_VIC01_48hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_VIC01_48hpi_child_data.filt_1 <- NormalizeData(Short_read_VIC01_48hpi_child_data.filt_1)
Short_read_VIC01_48hpi_child_data.filt_2 <- NormalizeData(Short_read_VIC01_48hpi_child_data.filt_2)
Short_read_VIC01_48hpi_child_data.filt_3 <- NormalizeData(Short_read_VIC01_48hpi_child_data.filt_3)



Short_read_VIC01_48hpi_child_data.filt_1 <- FindVariableFeatures(Short_read_VIC01_48hpi_child_data.filt_1, selection.method = "vst")
Short_read_VIC01_48hpi_child_data.filt_2 <- FindVariableFeatures(Short_read_VIC01_48hpi_child_data.filt_2, selection.method = "vst")
Short_read_VIC01_48hpi_child_data.filt_3 <- FindVariableFeatures(Short_read_VIC01_48hpi_child_data.filt_3, selection.method = "vst")


Short_read_VIC01_48hpi_child_data.filt_1 <- ScaleData(Short_read_VIC01_48hpi_child_data.filt_1, features = rownames(Short_read_VIC01_48hpi_child_data.filt_1))
Short_read_VIC01_48hpi_child_data.filt_2 <- ScaleData(Short_read_VIC01_48hpi_child_data.filt_2, features = rownames(Short_read_VIC01_48hpi_child_data.filt_2))
Short_read_VIC01_48hpi_child_data.filt_3 <- ScaleData(Short_read_VIC01_48hpi_child_data.filt_3, features = rownames(Short_read_VIC01_48hpi_child_data.filt_3))



Short_read_VIC01_48hpi_child_data.filt_1 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_1, features = VariableFeatures(Short_read_VIC01_48hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_48hpi_child_data.filt_2 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_2, features = VariableFeatures(Short_read_VIC01_48hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_48hpi_child_data.filt_3 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_3, features = VariableFeatures(Short_read_VIC01_48hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)



Short_read_VIC01_48hpi_child_data.filt_1 <- CellCycleScoring(Short_read_VIC01_48hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_48hpi_child_data.filt_2 <- CellCycleScoring(Short_read_VIC01_48hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_48hpi_child_data.filt_3 <- CellCycleScoring(Short_read_VIC01_48hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)




#Alternate Workflow regression
Short_read_VIC01_48hpi_child_data.filt_1$CC.Difference <- Short_read_VIC01_48hpi_child_data.filt_1$S.Score - Short_read_VIC01_48hpi_child_data.filt_1$G2M.Score
Short_read_VIC01_48hpi_child_data.filt_2$CC.Difference <- Short_read_VIC01_48hpi_child_data.filt_2$S.Score - Short_read_VIC01_48hpi_child_data.filt_2$G2M.Score
Short_read_VIC01_48hpi_child_data.filt_3$CC.Difference <- Short_read_VIC01_48hpi_child_data.filt_3$S.Score - Short_read_VIC01_48hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_VIC01_48hpi_child_data.filt_1 <- SCTransform(Short_read_VIC01_48hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_48hpi_child_data.filt_2 <- SCTransform(Short_read_VIC01_48hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_48hpi_child_data.filt_3 <- SCTransform(Short_read_VIC01_48hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)

Short_read_VIC01_48hpi_child_data.filt_1 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_48hpi_child_data.filt_2 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_48hpi_child_data.filt_3 <- RunPCA(Short_read_VIC01_48hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)

DimPlot(Short_read_VIC01_48hpi_child_data.filt_1)
DimPlot(Short_read_VIC01_48hpi_child_data.filt_2)
DimPlot(Short_read_VIC01_48hpi_child_data.filt_3)


# store mitochondrial percentage in object meta data
Short_read_VIC01_72hpi_child<- PercentageFeatureSet(Short_read_VIC01_72hpi_child, pattern = "^MT-", col.name = "percent_mito")
#Short_read_VIC01_72hpi_child_f <- PercentageFeatureSet(Short_read_VIC01_72hpi_child_f, pattern = "^MT-", col.name = "percent_mito")

# store ribosomal percentage in object meta data
Short_read_VIC01_72hpi_child<- PercentageFeatureSet(Short_read_VIC01_72hpi_child, "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
Short_read_VIC01_72hpi_child<- PercentageFeatureSet(Short_read_VIC01_72hpi_child, "^HB[^(P)]", col.name = "percent_hb")


# run sctransform (with glmGamPoi)
#Short_read_VIC01_72hpi_child<- SCTransform(Short_read_VIC01_72hpi_child, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)
#Short_read_VIC01_72hpi_child_f <- SCTransform(Short_read_VIC01_72hpi_child_f, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = FALSE)

#Filter only Singlet

Short_read_VIC01_72hpi_child_data <- subset(Short_read_VIC01_72hpi_child, subset = demux.doublet.call =="SNG")


#detection based testing
Short_read_VIC01_72hpi_child_selected_c <- WhichCells(Short_read_VIC01_72hpi_child_data, expression = nFeature_RNA > 200 & nFeature_RNA < 9000)
Short_read_VIC01_72hpi_child_selected_f <- rownames(Short_read_VIC01_72hpi_child_data)[Matrix::rowSums(Short_read_VIC01_72hpi_child_data) > 3]

Short_read_VIC01_72hpi_child_data.filt <- subset(Short_read_VIC01_72hpi_child_data, features = Short_read_VIC01_72hpi_child_selected_f, cells = Short_read_VIC01_72hpi_child_selected_c)
dim(Short_read_VIC01_72hpi_child_data.filt)




#mito/ribo filtering
Short_read_VIC01_72hpi_child_selected_mito <- WhichCells(Short_read_VIC01_72hpi_child_data.filt, expression = percent_mito < 20)
Short_read_VIC01_72hpi_child_selected_ribo <- WhichCells(Short_read_VIC01_72hpi_child_data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
Short_read_VIC01_72hpi_child_data.filt <- subset(Short_read_VIC01_72hpi_child_data.filt, cells = Short_read_VIC01_72hpi_child_selected_mito)
Short_read_VIC01_72hpi_child_data.filt <- subset(Short_read_VIC01_72hpi_child_data.filt, cells = Short_read_VIC01_72hpi_child_selected_ribo)

dim(Short_read_VIC01_72hpi_child_data.filt)

table(Short_read_VIC01_72hpi_child_data.filt$orig.ident)


Short_read_VIC01_72hpi_child_data.filt_1 <- subset(Short_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample4_Sample4")
Short_read_VIC01_72hpi_child_data.filt_2 <- subset(Short_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample5_Sample5")
Short_read_VIC01_72hpi_child_data.filt_3 <- subset(Short_read_VIC01_72hpi_child_data.filt, subset = Sample =="Sample6_Sample6")


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


Short_read_VIC01_72hpi_child_data.filt_1 <- NormalizeData(Short_read_VIC01_72hpi_child_data.filt_1)
Short_read_VIC01_72hpi_child_data.filt_2 <- NormalizeData(Short_read_VIC01_72hpi_child_data.filt_2)
Short_read_VIC01_72hpi_child_data.filt_3 <- NormalizeData(Short_read_VIC01_72hpi_child_data.filt_3)



Short_read_VIC01_72hpi_child_data.filt_1 <- FindVariableFeatures(Short_read_VIC01_72hpi_child_data.filt_1, selection.method = "vst")
Short_read_VIC01_72hpi_child_data.filt_2 <- FindVariableFeatures(Short_read_VIC01_72hpi_child_data.filt_2, selection.method = "vst")
Short_read_VIC01_72hpi_child_data.filt_3 <- FindVariableFeatures(Short_read_VIC01_72hpi_child_data.filt_3, selection.method = "vst")


Short_read_VIC01_72hpi_child_data.filt_1 <- ScaleData(Short_read_VIC01_72hpi_child_data.filt_1, features = rownames(Short_read_VIC01_72hpi_child_data.filt_1))
Short_read_VIC01_72hpi_child_data.filt_2 <- ScaleData(Short_read_VIC01_72hpi_child_data.filt_2, features = rownames(Short_read_VIC01_72hpi_child_data.filt_2))
Short_read_VIC01_72hpi_child_data.filt_3 <- ScaleData(Short_read_VIC01_72hpi_child_data.filt_3, features = rownames(Short_read_VIC01_72hpi_child_data.filt_3))



Short_read_VIC01_72hpi_child_data.filt_1 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_1, features = VariableFeatures(Short_read_VIC01_72hpi_child_data.filt_1), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_72hpi_child_data.filt_2 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_2, features = VariableFeatures(Short_read_VIC01_72hpi_child_data.filt_2), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
Short_read_VIC01_72hpi_child_data.filt_3 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_3, features = VariableFeatures(Short_read_VIC01_72hpi_child_data.filt_3), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)


Short_read_VIC01_72hpi_child_data.filt_1 <- CellCycleScoring(Short_read_VIC01_72hpi_child_data.filt_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_72hpi_child_data.filt_2 <- CellCycleScoring(Short_read_VIC01_72hpi_child_data.filt_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Short_read_VIC01_72hpi_child_data.filt_3 <- CellCycleScoring(Short_read_VIC01_72hpi_child_data.filt_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)





#Alternate Workflow regression
Short_read_VIC01_72hpi_child_data.filt_1$CC.Difference <- Short_read_VIC01_72hpi_child_data.filt_1$S.Score - Short_read_VIC01_72hpi_child_data.filt_1$G2M.Score
Short_read_VIC01_72hpi_child_data.filt_2$CC.Difference <- Short_read_VIC01_72hpi_child_data.filt_2$S.Score - Short_read_VIC01_72hpi_child_data.filt_2$G2M.Score
Short_read_VIC01_72hpi_child_data.filt_3$CC.Difference <- Short_read_VIC01_72hpi_child_data.filt_3$S.Score - Short_read_VIC01_72hpi_child_data.filt_3$G2M.Score

# run sctransform (with glmGamPoi)
Short_read_VIC01_72hpi_child_data.filt_1 <- SCTransform(Short_read_VIC01_72hpi_child_data.filt_1, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_72hpi_child_data.filt_2 <- SCTransform(Short_read_VIC01_72hpi_child_data.filt_2, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
Short_read_VIC01_72hpi_child_data.filt_3 <- SCTransform(Short_read_VIC01_72hpi_child_data.filt_3, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)


Short_read_VIC01_72hpi_child_data.filt_1 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_1, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_72hpi_child_data.filt_2 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_2, features = c(s.genes, g2m.genes),approx=FALSE)
Short_read_VIC01_72hpi_child_data.filt_3 <- RunPCA(Short_read_VIC01_72hpi_child_data.filt_3, features = c(s.genes, g2m.genes),approx=FALSE)


DimPlot(Short_read_VIC01_72hpi_child_data.filt_1)
DimPlot(Short_read_VIC01_72hpi_child_data.filt_2)
DimPlot(Short_read_VIC01_72hpi_child_data.filt_3)



save(Long_read_UK_72hpi_adult_data.filt_1, file = "Long_read_UK_72hpi_adult_data.filt_1.RData")
save(Long_read_uninfected_adult_data.filt_1, file =  "Long_read_uninfected_adult_data.filt_1.RData")
save(Long_read_VIC01_48hpi_adult_data.filt_1, file = "Long_read_VIC01_48hpi_adult_data.filt_1.RData")
save(Long_read_VIC01_72hpi_adult_data.filt_1, file = "Long_read_VIC01_72hpi_adult_data.filt_1.RData")
save(Short_read_UK_72hpi_adult_data.filt_1, file = "Short_read_UK_72hpi_adult_data.filt_1.RData")
save(Short_read_uninfected_adult_data.filt_1, file = "Short_read_uninfected_adult_data.filt_1.RData")
save(Short_read_VIC01_48hpi_adult_data.filt_1, file = "Short_read_VIC01_48hpi_adult_data.filt_1.RData")
save(Short_read_VIC01_72hpi_adult_data.filt_1, file = "Short_read_VIC01_72hpi_adult_data.filt_1.RData")
save(Long_read_UK_72hpi_child_data.filt_1, file = "Long_read_UK_72hpi_child_data.filt_1.RData")
save(Long_read_uninfected_child_data.filt_1, file = "Long_read_uninfected_child_data.filt_1.RData")
save(Long_read_VIC01_48hpi_child_data.filt_1, file =  "Long_read_VIC01_48hpi_child_data.filt_1.RData")
save(Long_read_VIC01_72hpi_child_data.filt_1, file = "Long_read_VIC01_72hpi_child_data.filt_1.RData")
save(Short_read_UK_72hpi_child_data.filt_1, file = "Short_read_UK_72hpi_child_data.filt_1.RData")
save(Short_read_uninfected_child_data.filt_1, file = "Short_read_uninfected_child_data.filt_1.RData")
save(Short_read_VIC01_48hpi_child_data.filt_1, file = "Short_read_VIC01_48hpi_child_data.filt_1.RData")
save(Short_read_VIC01_72hpi_child_data.filt_1, file = "Short_read_VIC01_72hpi_child_data.filt_1.RData")


save(Long_read_UK_72hpi_adult_data.filt_2, file = "Long_read_UK_72hpi_adult_data.filt_2.RData")
save(Long_read_uninfected_adult_data.filt_2, file =  "Long_read_uninfected_adult_data.filt_2.RData")
save(Long_read_VIC01_48hpi_adult_data.filt_2, file = "Long_read_VIC01_48hpi_adult_data.filt_2.RData")
save(Long_read_VIC01_72hpi_adult_data.filt_2, file = "Long_read_VIC01_72hpi_adult_data.filt_2.RData")
save(Short_read_UK_72hpi_adult_data.filt_2, file = "Short_read_UK_72hpi_adult_data.filt_2.RData")
save(Short_read_uninfected_adult_data.filt_2, file = "Short_read_uninfected_adult_data.filt_2.RData")
save(Short_read_VIC01_48hpi_adult_data.filt_2, file = "Short_read_VIC01_48hpi_adult_data.filt_2.RData")
save(Short_read_VIC01_72hpi_adult_data.filt_2, file = "Short_read_VIC01_72hpi_adult_data.filt_2.RData")
save(Long_read_UK_72hpi_child_data.filt_2, file = "Long_read_UK_72hpi_child_data.filt_2.RData")
save(Long_read_uninfected_child_data.filt_2, file = "Long_read_uninfected_child_data.filt_2.RData")
save(Long_read_VIC01_48hpi_child_data.filt_2, file =  "Long_read_VIC01_48hpi_child_data.filt_2.RData")
save(Long_read_VIC01_72hpi_child_data.filt_2, file = "Long_read_VIC01_72hpi_child_data.filt_2.RData")
save(Short_read_UK_72hpi_child_data.filt_2, file = "Short_read_UK_72hpi_child_data.filt_2.RData")
save(Short_read_uninfected_child_data.filt_2, file = "Short_read_uninfected_child_data.filt_2.RData")
save(Short_read_VIC01_48hpi_child_data.filt_2, file = "Short_read_VIC01_48hpi_child_data.filt_2.RData")
save(Short_read_VIC01_72hpi_child_data.filt_2, file = "Short_read_VIC01_72hpi_child_data.filt_2.RData")


save(Long_read_UK_72hpi_adult_data.filt_3, file = "Long_read_UK_72hpi_adult_data.filt_3.RData")
save(Long_read_uninfected_adult_data.filt_3, file =  "Long_read_uninfected_adult_data.filt_3.RData")
save(Long_read_VIC01_48hpi_adult_data.filt_3, file = "Long_read_VIC01_48hpi_adult_data.filt_3.RData")
save(Long_read_VIC01_72hpi_adult_data.filt_3, file = "Long_read_VIC01_72hpi_adult_data.filt_3.RData")
save(Short_read_UK_72hpi_adult_data.filt_3, file = "Short_read_UK_72hpi_adult_data.filt_3.RData")
save(Short_read_uninfected_adult_data.filt_3, file = "Short_read_uninfected_adult_data.filt_3.RData")
save(Short_read_VIC01_48hpi_adult_data.filt_3, file = "Short_read_VIC01_48hpi_adult_data.filt_3.RData")
save(Short_read_VIC01_72hpi_adult_data.filt_3, file = "Short_read_VIC01_72hpi_adult_data.filt_3.RData")
save(Long_read_UK_72hpi_child_data.filt_3, file = "Long_read_UK_72hpi_child_data.filt_3.RData")
save(Long_read_uninfected_child_data.filt_3, file = "Long_read_uninfected_child_data.filt_3.RData")
save(Long_read_VIC01_48hpi_child_data.filt_3, file =  "Long_read_VIC01_48hpi_child_data.filt_3.RData")
save(Long_read_VIC01_72hpi_child_data.filt_3, file = "Long_read_VIC01_72hpi_child_data.filt_3.RData")
save(Short_read_UK_72hpi_child_data.filt_3, file = "Short_read_UK_72hpi_child_data.filt_3.RData")
save(Short_read_uninfected_child_data.filt_3, file = "Short_read_uninfected_child_data.filt_3.RData")
save(Short_read_VIC01_48hpi_child_data.filt_3, file = "Short_read_VIC01_48hpi_child_data.filt_3.RData")
save(Short_read_VIC01_72hpi_child_data.filt_3, file = "Short_read_VIC01_72hpi_child_data.filt_3.RData")

