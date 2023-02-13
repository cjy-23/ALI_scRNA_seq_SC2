#setwd("Z:/User/Jessie/Projects/scRNA-seq/Analysis/Seurat/seurat_qc_3")
#setwd("C:/Users/jessiejieyou/Documents/scRNA_seq")

library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(metap)
library(multtest)
library(dplyr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(RCurl)
library(cowplot)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(dittoSeq)
#memory.limit(size=1600000)

load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_adult_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_child_data.filt_1.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_child_data.filt_1.RData")


load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_adult_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_child_data.filt_2.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_child_data.filt_2.RData")


load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_adult_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_UK_72hpi_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_uninfected_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_48hpi_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Long_read_VIC01_72hpi_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_UK_72hpi_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_uninfected_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_48hpi_child_data.filt_3.RData")
load("/data/gpfs/projects/punim1466/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_1_2_adult_child_inf/Short_read_VIC01_72hpi_child_data.filt_3.RData")




##---------merge
immune.combined.sct <- merge(Long_read_UK_72hpi_adult_data.filt_1, c(Long_read_UK_72hpi_adult_data.filt_2,Long_read_UK_72hpi_adult_data.filt_3,
Long_read_uninfected_adult_data.filt_1,Long_read_uninfected_adult_data.filt_2,Long_read_uninfected_adult_data.filt_3,
Long_read_VIC01_48hpi_adult_data.filt_1,Long_read_VIC01_48hpi_adult_data.filt_2,Long_read_VIC01_48hpi_adult_data.filt_3,
Long_read_VIC01_72hpi_adult_data.filt_1, Long_read_VIC01_72hpi_adult_data.filt_2, Long_read_VIC01_72hpi_adult_data.filt_3,
Short_read_UK_72hpi_adult_data.filt_1,  Short_read_UK_72hpi_adult_data.filt_2,  Short_read_UK_72hpi_adult_data.filt_3,
Short_read_uninfected_adult_data.filt_1,  Short_read_uninfected_adult_data.filt_2,  Short_read_uninfected_adult_data.filt_3,
Short_read_VIC01_48hpi_adult_data.filt_1,  Short_read_VIC01_48hpi_adult_data.filt_2,  Short_read_VIC01_48hpi_adult_data.filt_3,
Short_read_VIC01_72hpi_adult_data.filt_1,  Short_read_VIC01_72hpi_adult_data.filt_2,  Short_read_VIC01_72hpi_adult_data.filt_3,
Long_read_UK_72hpi_child_data.filt_1, Long_read_UK_72hpi_child_data.filt_2,Long_read_UK_72hpi_child_data.filt_3,
Long_read_uninfected_child_data.filt_1,Long_read_uninfected_child_data.filt_2,Long_read_uninfected_child_data.filt_3,
Long_read_VIC01_48hpi_child_data.filt_1,Long_read_VIC01_48hpi_child_data.filt_2,Long_read_VIC01_48hpi_child_data.filt_3,
Long_read_VIC01_72hpi_child_data.filt_1, Long_read_VIC01_72hpi_child_data.filt_2, Long_read_VIC01_72hpi_child_data.filt_3,
Short_read_UK_72hpi_child_data.filt_1,  Short_read_UK_72hpi_child_data.filt_2,  Short_read_UK_72hpi_child_data.filt_3,
Short_read_uninfected_child_data.filt_1,  Short_read_uninfected_child_data.filt_2,  Short_read_uninfected_child_data.filt_3,
Short_read_VIC01_48hpi_child_data.filt_1,  Short_read_VIC01_48hpi_child_data.filt_2,  Short_read_VIC01_48hpi_child_data.filt_3,
Short_read_VIC01_72hpi_child_data.filt_1,  Short_read_VIC01_72hpi_child_data.filt_2,  Short_read_VIC01_72hpi_child_data.filt_3))


##https://github.com/satijalab/seurat/issues/2814

VariableFeatures(immune.combined.sct[["SCT"]]) <- rownames(immune.combined.sct[["SCT"]]@scale.data)

save(immune.combined.sct, file = "immune.combined.sct.RData")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.5)


# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "type", raster=FALSE)
pdf(file="type.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p1
dev.off()

p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = FALSE, raster=FALSE)
pdf(file="cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p2
dev.off()


pdf(file="type_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p1+p2
dev.off()


p3 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "group", raster=FALSE)
pdf(file="group.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p3
dev.off()

pdf(file="group_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p3+p2
dev.off()


p5 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Sample", raster=FALSE)
pdf(file="donor.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p5
dev.off()

pdf(file="donor_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p5+p2
dev.off()

p6 <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "Sample", raster=FALSE)
pdf(file="donor_split.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p6
dev.off()

p7 <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "group", raster=FALSE)
pdf(file="group_split.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p7
dev.off()



p8 <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "type", raster=FALSE)
pdf(file="type_split.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p8
dev.off()

p9 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Infection", raster=FALSE)
pdf(file="Infection.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p9
dev.off()

pdf(file="infection_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)

p9+p2
dev.off()


p10 <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "Infection", raster=FALSE)
pdf(file="infection_split.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p10
dev.off()


p11 <- FeaturePlot(immune.combined.sct, reduction = "umap", features = "percent_mito", raster=FALSE)
pdf(file="mito.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p11
dev.off()

pdf(file="mito_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p11 + p2
dev.off()


p12 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Phase")
pdf(file="ccphase.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p12
dev.off()


pdf(file="phase_cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p12+p2
dev.off()