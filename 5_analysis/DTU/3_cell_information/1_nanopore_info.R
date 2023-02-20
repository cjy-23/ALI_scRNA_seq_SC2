library(Seurat)
library(dplyr)
library(purrr)
library(tibble)
library(dplyr)
library(Libra)
library(xlsx)
library(stringr)

setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/")
#setwd("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/seurat_meta_data_2/")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_v2/v1_notidy/split_1_4_9_10_12/sub_RNA_5.RData")
#scp_meta <- read.delim("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_2_inf/dim_20/reso_0.3_dim_20/cluster_7/metadata_donor_6_no_int.tsv")

#-------------Assign names to clusters
Idents(sub_5) <- sub_5$test_5
##redo thissss
sub_5 <- RenameIdents(sub_5, 
                      
                      '0'='Suprabasal',
                      '1_0'='Secretory-1',
                      '1_1'='Secretory-1',
                      '1_2'='Secretory-1',
                      '1_3'='Secretory-1',
                      '1_4'='Secretory-1',
                      '1_5'='Goblet',
                      '1_6'='Secretory-1',
                      '1_7'='Goblet-Ionocyte',
                      '2'	="Basal-1",
                      '3'	="Ciliated-1",
                      '4_0'	="Goblet",
                      '4_1'	="Secretory-2",
                      '4_2'	="Goblet",
                      '4_3'	="Secretory-2",
                      '4_4'	="Secretory-2",
                      '4_5'	="Secretory-2",
                      '4_6'	="Secretory-2",
                      '4_7'	="Secretory-2",
                      '4_8'	="Secretory-2",
                      '5'	="Ciliated-3",
                      '6'	="Cycling-basal",
                      '7'	="Ciliated-2",
                      '8'	="Basal-2",
                      '9_0'	="Ciliated-SC2",
                      '9_1'	="Secretory-Ciliated-SC2",
                      '10_0' ="Goblet-IFN-stim",
                      '10_1' ="Secretory-IFN-stim",
                      '10_2' ="Secretory-IFN-stim",
                      '10_3' ="Secretory-IFN-stim",
                      '10_4' ="Goblet-IFN-stim",
                      '10_5' ="Secretory-IFN-stim",
                      '10_6' ="Secretory-IFN-stim",
                      '10_7' ="Secretory-IFN-stim",
                      '11'	="Deuterosomal",
                      '12_0' = "Ionocyte",
                      '12_1' = "Brush/Tuft",
                      '13' = "Secretory-Ciliated"
)
sub_5$celltype <- Idents(sub_5)
data <- sub_5@meta.data

#cell barcode add
data$barcode <- rownames(data)

#data remove barcode suffix
data$barcode <- str_sub(data$barcode,1,nchar(data$barcode)-4)


#Extract cell barcode
#barcode <- rownames(data)
#barcode <- str_sub(rname,1,nchar(barcode)-4)
#barcode <- as.data.frame.character(barcode)
#d <- duplicated(rname)
#d <- as.data.frame(d)
#Duplicated cell-barcodes
#d_2 <- subset(d, d=="TRUE")


#Add barcode data to meta.data
#data$barcode <- merge (data,barcode)
#data$barcode <- as.character(data$barcode)

#Subset Long_read_UK_72hpi_adult

Long_read_UK_72hpi_adult <- subset(data, type=="Long_read_UK_72hpi_adult")
u <- duplicated(Long_read_UK_72hpi_adult$barcode)
u <- as.data.frame(u) #No duplications

Long_read_uninfected_adult <- subset(data, type=="Long_read_uninfected_adult")
p <- duplicated(Long_read_uninfected_adult$barcode)
p <- as.data.frame(p) #No duplications

Long_read_UK_72hpi_child <- subset(data, type=="Long_read_UK_72hpi_child")
l <- duplicated(Long_read_UK_72hpi_child$barcode)
l <- as.data.frame(l) #No duplications

Long_read_uninfected_child <- subset(data, type=="Long_read_uninfected_child")
i <- duplicated(Long_read_uninfected_child$barcode)
i <- as.data.frame(i) #No duplications

Long_read_VIC01_72hpi_adult <- subset(data, type=="Long_read_VIC01_72hpi_adult")
c <- duplicated(Long_read_VIC01_72hpi_adult$barcode)
c <- as.data.frame(c) #No duplications

Long_read_VIC01_72hpi_child <- subset(data, type=="Long_read_VIC01_72hpi_child")
a <- duplicated(Long_read_VIC01_72hpi_child$barcode)
a <- as.data.frame(a) #No duplications


Long_read_VIC01_48hpi_adult <- subset(data, type=="Long_read_VIC01_48hpi_adult")
t <- duplicated(Long_read_VIC01_48hpi_adult$barcode)
t <- as.data.frame(t) #No duplications

Long_read_VIC01_48hpi_child <- subset(data, type=="Long_read_VIC01_48hpi_child")
e <- duplicated(Long_read_VIC01_48hpi_child$barcode)
e <- as.data.frame(e) #No duplications



#Save data
write.csv(Long_read_UK_72hpi_adult, file="Long_read_UK_72hpi_adult_barcode.csv", row.names=TRUE)
write.csv(Long_read_uninfected_adult, file="Long_read_uninfected_adult_barcode.csv", row.names=TRUE)
write.csv(Long_read_UK_72hpi_child, file="Long_read_UK_72hpi_child_barcode.csv", row.names=TRUE)
write.csv(Long_read_uninfected_child, file="Long_read_uninfected_child_barcode.csv", row.names=TRUE)
write.csv(Long_read_VIC01_72hpi_adult, file="Long_read_VIC01_72hpi_adult_barcode.csv", row.names=TRUE)
write.csv(Long_read_VIC01_72hpi_child, file="Long_read_VIC01_72hpi_child_barcode.csv", row.names=TRUE)
write.csv(Long_read_VIC01_48hpi_adult, file="Long_read_VIC01_48hpi_adult_barcode.csv", row.names=TRUE)
write.csv(Long_read_VIC01_48hpi_child, file="Long_read_VIC01_48hpi_child_barcode.csv", row.names=TRUE)

#Subset data just by Secretory_Ciliated_1 cells
#Long_read_UK_72hpi_adult_sc1 <- subset(Long_read_UK_72hpi_adult, celltype=="Secretory_Ciliated_1") #276 cells
#Long_read_uninfected_adult_sc1 <- subset(Long_read_uninfected_adult, celltype=="Secretory_Ciliated_1") #100 cells

#Extract barcode
#Long_read_UK_72hpi_adult_sc1_bc <- Long_read_UK_72hpi_adult_sc1$barcode
#Long_read_uninfected_adult_sc1_bc <- Long_read_uninfected_adult_sc1$barcode

#Save data
#write.csv(Long_read_UK_72hpi_adult_sc1_bc, file="Long_read_UK_72hpi_adult_sc1_bc.csv")
#write.csv(Long_read_uninfected_adult_sc1_bc, file="Long_read_uninfected_adult_sc1_bc.csv")




#sub_5$celltype.group <- paste(Idents(sub_5), sub_5$group, sub_5$Infection, sep = "_")
#sub_5$celltype <- Idents(sub_5)
#Idents(sub_5) <- "celltype.group"


#------------------------------ Merge data from seurat meta data and scp_meta.txt
libA <- Long_read_UK_72hpi_adult[, c(34, 33, 4, 15, 7)]
libC <- Long_read_uninfected_adult[, c(34, 33, 4, 15, 7)]
libB <- Long_read_UK_72hpi_child[, c(34, 33, 4, 15, 7)]
libD <- Long_read_uninfected_child[, c(34, 33, 4, 15, 7)]
libE <- Long_read_VIC01_48hpi_adult[, c(34, 33, 4, 15, 7)]
libF <- Long_read_VIC01_48hpi_child[, c(34, 33, 4, 15, 7)]
libG <- Long_read_VIC01_72hpi_adult[, c(34, 33, 4, 15, 7)]
libH <- Long_read_VIC01_72hpi_child[, c(34, 33, 4, 15, 7)]




colnames(libA) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libB) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libC) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libD) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libE) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libF) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libG) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")
colnames(libH) <- c("NAME", "Cell_Type", "Cell_State", "Cohort", "Patient")

#Remove rownames
rownames(libA)<-NULL
rownames(libB)<-NULL
rownames(libC)<-NULL
rownames(libD)<-NULL
rownames(libE)<-NULL
rownames(libF)<-NULL
rownames(libG)<-NULL
rownames(libH)<-NULL


libB$NAME <- str_sub(libB$NAME,1,nchar(libB$NAME)-1)
libD$NAME <- str_sub(libD$NAME,1,nchar(libD$NAME)-1)
libF$NAME <- str_sub(libF$NAME,1,nchar(libF$NAME)-1)
libG$NAME <- str_sub(libG$NAME,1,nchar(libG$NAME)-1)
libH$NAME <- str_sub(libH$NAME,1,nchar(libH$NAME)-1)




write.table(libA, file="libA_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libB, file="libB_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libC, file="libC_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libD, file="libD_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libE, file="libE_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libF, file="libF_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libG, file="libG_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(libH, file="libH_scp_meta.txt", sep="\t", quote=FALSE, row.names=FALSE)
