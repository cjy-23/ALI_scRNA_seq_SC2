## Load packages
library(loomR)
library(hdf5r)
library(R6)
library(iterators)
library(itertools)
library(rhdf5)
library(dplyr)
library(SeuratWrappers)
library(pcaMethods)
library(Seurat)
#BiocManager::install("pcaMethods")



setwd("/velocity/velocyto_filtered/")

#1. Copy the files from scratch, upload to mediaflux, 

##Merge
Long_read_UK_72hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Long_read_UK_72hpi_adult/velocyto/possorted_genome_bam_G0BQ9.loom")
Long_read_UK_72hpi_child <- ReadVelocity(file = "/velocyto_filtered/Long_read_UK_72hpi_child/velocyto/possorted_genome_bam_VI1LA.loom")
Long_read_uninfected_adult <- ReadVelocity(file = "/velocyto_filtered/Long_read_uninfected_adult/velocyto/possorted_genome_bam_3913P.loom")
Long_read_uninfected_child <- ReadVelocity(file = "/velocyto_filtered/Long_read_uninfected_child/velocyto/possorted_genome_bam_OC15V.loom")
Long_read_VIC01_48hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Long_read_VIC01_48hpi_adult/velocyto/possorted_genome_bam_EYWGZ.loom")
Long_read_VIC01_48hpi_child <- ReadVelocity(file = "/velocyto_filtered/Long_read_VIC01_48hpi_child/velocyto/possorted_genome_bam_U1DUE.loom")
Long_read_VIC01_72hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Long_read_VIC01_72hpi_adult/velocyto/possorted_genome_bam_OYCZC.loom")
Long_read_VIC01_72hpi_child <- ReadVelocity(file = "/velocyto_filtered/Long_read_VIC01_72hpi_child/velocyto/possorted_genome_bam_A1C6D.loom")

                                             
Short_read_UK_72hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Short_read_UK_72hpi_adult/velocyto/possorted_genome_bam_TASCN.loom")
Short_read_UK_72hpi_child <- ReadVelocity(file = "/velocyto_filtered/Short_read_UK_72hpi_child/velocyto/possorted_genome_bam_G5UGR.loom")
Short_read_uninfected_adult <- ReadVelocity(file = "/velocyto_filtered/Short_read_uninfected_adult/velocyto/possorted_genome_bam_84BX7.loom")
Short_read_uninfected_child <- ReadVelocity(file = "/velocyto_filtered/Short_read_uninfected_child/velocyto/possorted_genome_bam_B7UTP.loom")
Short_read_VIC01_48hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Short_read_VIC01_48hpi_adult/velocyto/possorted_genome_bam_8MSW8.loom")
Short_read_VIC01_48hpi_child <- ReadVelocity(file = "/velocyto_filtered/Short_read_VIC01_48hpi_child/velocyto/possorted_genome_bam_8F3WC.loom")
Short_read_VIC01_72hpi_adult <- ReadVelocity(file = "/velocyto_filtered/Short_read_VIC01_72hpi_adult/velocyto/possorted_genome_bam_0AK2O.loom")
Short_read_VIC01_72hpi_child <- ReadVelocity(file = "/velocyto_filtered/Short_read_VIC01_72hpi_child/velocyto/possorted_genome_bam_2HAMB.loom")




##Convert to Seurat object

Long_read_UK_72hpi_adult_seurat <- as.Seurat(x = Long_read_UK_72hpi_adult)
Long_read_UK_72hpi_child_seurat <- as.Seurat(x = Long_read_UK_72hpi_child)
Long_read_uninfected_adult_seurat <- as.Seurat(x = Long_read_uninfected_adult)
Long_read_uninfected_child_seurat <- as.Seurat(x = Long_read_uninfected_child)
Long_read_VIC01_48hpi_adult_seurat <- as.Seurat(x = Long_read_VIC01_48hpi_adult)
Long_read_VIC01_48hpi_child_seurat <- as.Seurat(x = Long_read_VIC01_48hpi_child)
Long_read_VIC01_72hpi_adult_seurat <- as.Seurat(x = Long_read_VIC01_72hpi_adult)
Long_read_VIC01_72hpi_child_seurat <- as.Seurat(x = Long_read_VIC01_72hpi_child)


Short_read_UK_72hpi_adult_seurat <- as.Seurat(x = Short_read_UK_72hpi_adult)
Short_read_UK_72hpi_child_seurat <- as.Seurat(x = Short_read_UK_72hpi_child)
Short_read_uninfected_adult_seurat <- as.Seurat(x = Short_read_uninfected_adult)
Short_read_uninfected_child_seurat <- as.Seurat(x = Short_read_uninfected_child)
Short_read_VIC01_48hpi_adult_seurat <- as.Seurat(x = Short_read_VIC01_48hpi_adult)
Short_read_VIC01_48hpi_child_seurat <- as.Seurat(x = Short_read_VIC01_48hpi_child)
Short_read_VIC01_72hpi_adult_seurat <- as.Seurat(x = Short_read_VIC01_72hpi_adult)
Short_read_VIC01_72hpi_child_seurat <- as.Seurat(x = Short_read_VIC01_72hpi_child)


#Add Type


Long_read_UK_72hpi_adult_seurat@meta.data$Type <- "Long_read_UK_72hpi_adult"
Long_read_UK_72hpi_child_seurat@meta.data$Type <- "Long_read_UK_72hpi_child"
Long_read_uninfected_adult_seurat@meta.data$Type <- "Long_read_uninfected_adult"
Long_read_uninfected_child_seurat@meta.data$Type <- "Long_read_uninfected_child"
Long_read_VIC01_48hpi_adult_seurat@meta.data$Type <- "Long_read_VIC01_48hpi_adult"
Long_read_VIC01_48hpi_child_seurat@meta.data$Type <- "Long_read_VIC01_48hpi_child"
Long_read_VIC01_72hpi_adult_seurat@meta.data$Type <- "Long_read_VIC01_72hpi_adult"
Long_read_VIC01_72hpi_child_seurat@meta.data$Type <- "Long_read_VIC01_72hpi_child"


Short_read_UK_72hpi_adult_seurat@meta.data$Type  <- "Short_read_UK_72hpi_adult"
Short_read_UK_72hpi_child_seurat@meta.data$Type  <- "Short_read_UK_72hpi_child"
Short_read_uninfected_adult_seurat@meta.data$Type  <- "Short_read_uninfected_adult"
Short_read_uninfected_child_seurat@meta.data$Type  <- "Short_read_uninfected_child"
Short_read_VIC01_48hpi_adult_seurat@meta.data$Type  <- "Short_read_VIC01_48hpi_adult"
Short_read_VIC01_48hpi_child_seurat@meta.data$Type  <- "Short_read_VIC01_48hpi_child"
Short_read_VIC01_72hpi_adult_seurat@meta.data$Type  <- "Short_read_VIC01_72hpi_adult"
Short_read_VIC01_72hpi_child_seurat@meta.data$Type  <- "Short_read_VIC01_72hpi_child"


#Add Group


Long_read_UK_72hpi_adult_seurat@meta.data$Group <- "UK_72hpi_adult"
Long_read_UK_72hpi_child_seurat@meta.data$Group <- "UK_72hpi_child"
Long_read_uninfected_adult_seurat@meta.data$Group <- "uninfected_adult"
Long_read_uninfected_child_seurat@meta.data$Group <- "uninfected_child"
Long_read_VIC01_48hpi_adult_seurat@meta.data$Group <- "VIC01_48hpi_adult"
Long_read_VIC01_48hpi_child_seurat@meta.data$Group <- "VIC01_48hpi_child"
Long_read_VIC01_72hpi_adult_seurat@meta.data$Group <- "VIC01_72hpi_adult"
Long_read_VIC01_72hpi_child_seurat@meta.data$Group <- "VIC01_72hpi_child"


Short_read_UK_72hpi_adult_seurat@meta.data$Group  <- "UK_72hpi_adult"
Short_read_UK_72hpi_child_seurat@meta.data$Group  <- "UK_72hpi_child"
Short_read_uninfected_adult_seurat@meta.data$Group  <- "uninfected_adult"
Short_read_uninfected_child_seurat@meta.data$Group  <- "uninfected_child"
Short_read_VIC01_48hpi_adult_seurat@meta.data$Group  <- "VIC01_48hpi_adult"
Short_read_VIC01_48hpi_child_seurat@meta.data$Group  <- "VIC01_48hpi_child"
Short_read_VIC01_72hpi_adult_seurat@meta.data$Group  <- "VIC01_72hpi_adult"
Short_read_VIC01_72hpi_child_seurat@meta.data$Group  <- "VIC01_72hpi_child"




#Test where are rows dup

feat <- as.data.frame(Long_read_UK_72hpi_adult[["spliced"]]@Dimnames[[1]])
feat2 <-as.data.frame(Long_read_UK_72hpi_adult[["spliced"]]@Dimnames[[2]])

dup1 <- which(duplicated(feat$`Long_read_UK_72hpi_adult[["spliced"]]@Dimnames[[1]]`))
which(duplicated(feat2$`Long_read_UK_72hpi_adult[["spliced"]]@Dimnames[[2]]`))
dups <- as.data.frame(feat[dup1,])




#make a loop for this-


  
list_1 <- list(Long_read_UK_72hpi_adult, Long_read_UK_72hpi_child, Long_read_uninfected_adult, Long_read_uninfected_child,
                                         Long_read_VIC01_48hpi_adult,Long_read_VIC01_48hpi_child, Long_read_VIC01_72hpi_adult,Long_read_VIC01_72hpi_child,
                                         Short_read_UK_72hpi_adult, Short_read_UK_72hpi_child, Short_read_uninfected_adult, Short_read_uninfected_child,
                                         Short_read_VIC01_48hpi_adult, Short_read_VIC01_48hpi_child, Short_read_VIC01_72hpi_adult, Short_read_VIC01_72hpi_child)



fun1 <- function(x) {
  feat <- as.data.frame(x[["spliced"]]@Dimnames[[1]])
  dup1 <- which(duplicated(feat$`x[["spliced"]]@Dimnames[[1]]`))
  dups <- as.data.frame(feat[dup1,])
  dups
  
}

list_2 <- lapply(list_1, fun1)


feat <- as.data.frame(Long_read_uninfected_adult[["spliced"]]@Dimnames[[1]])
dup1 <- which(duplicated(feat$`Long_read_uninfected_adult[["spliced"]]@Dimnames[[1]]`))
dups <- as.data.frame(feat[dup1,])


feat <- as.data.frame(Short_read_VIC01_48hpi_adult[["spliced"]]@Dimnames[[1]])
dup1 <- which(duplicated(feat$`Short_read_VIC01_48hpi_adult[["spliced"]]@Dimnames[[1]]`))
dups <- as.data.frame(feat[dup1,])




#find equivalent in seurat object
Long_read_UK_72hpi_adult_seurat
revised_genes <- as.data.frame(Long_read_UK_72hpi_adult_seurat@assays[["spliced"]]@counts@Dimnames[[1]])




## Merge data
merged_ob <- merge(x = Long_read_UK_72hpi_adult_seurat, y = c(Long_read_UK_72hpi_child_seurat, Long_read_uninfected_adult_seurat, Long_read_uninfected_child_seurat,
                                                              Long_read_VIC01_48hpi_adult_seurat,Long_read_VIC01_48hpi_child_seurat, Long_read_VIC01_72hpi_adult_seurat,Long_read_VIC01_72hpi_child_seurat,
                                                              Short_read_UK_72hpi_adult_seurat, Short_read_UK_72hpi_child_seurat, Short_read_uninfected_adult_seurat, Short_read_uninfected_child_seurat,
                                                              Short_read_VIC01_48hpi_adult_seurat, Short_read_VIC01_48hpi_child_seurat, Short_read_VIC01_72hpi_adult_seurat, Short_read_VIC01_72hpi_child_seurat), merge.data = TRUE)



##Change barcode name

merged_ob@meta.data$barcode_og <- rownames(merged_ob@meta.data)
merged_ob@meta.data$barcode_af <- gsub(".*:", "", merged_ob@meta.data$barcode_og)
merged_ob@meta.data$barcode_af1 <- gsub("x", "", merged_ob@meta.data$barcode_af)


#Save merged object

save(merged_ob, file="merged_subset.RData")
