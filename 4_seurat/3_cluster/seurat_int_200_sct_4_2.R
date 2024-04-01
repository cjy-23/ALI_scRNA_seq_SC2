library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
#library(metap)
#library(multtest)
#library(dplyr)
#library(data.table)
#library(stringr)
#library(dplyr)
#library(tidyverse)
#library(RCurl)
#library(cowplot)
#library(ggplot2)
#library(dplyr)
#library(hrbrthemes)
#library(viridis)
#library(dittoSeq)
#memory.limit(size=1600000)


load("/analysis/novaseq_virus_trim/seurat/adult_sample_split_sct_donor/sctransform_v2/v1_notidy/immune.combined.sct.RData")


immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE, approx=FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, dims = 1:20)
immune.combined.sct <- FindNeighbors(immune.combined.sct, dims = 1:20)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)



#immune.combined.sct@commands$FindClusters$graph.name


sub <- FindSubCluster(
        immune.combined.sct,
        12,
        graph.name="SCT_snn",
        subcluster.name = "test",
        resolution = 0.1,
        algorithm = 1
)

Idents(sub) <- "test"

sub_2 <- FindSubCluster(
        sub,
        "4",
        graph.name="SCT_snn",
        subcluster.name = "test_2",
        resolution = 0.5,
        algorithm = 1
)
Idents(sub_2) <- "test_2"



sub_3 <- FindSubCluster(
        sub_2,
        "9",
        graph.name="SCT_snn",
        subcluster.name = "test_3",
        resolution = 0.1,
        algorithm = 1
)
Idents(sub_3) <- "test_3"



sub_4 <- FindSubCluster(
        sub_3,
        "1",
        graph.name="SCT_snn",
        subcluster.name = "test_4",
        resolution = 0.4,
        algorithm = 1
)
Idents(sub_4) <- "test_4"


sub_5 <- FindSubCluster(
        sub_4,
        "10",
        graph.name="SCT_snn",
        subcluster.name = "test_5",
        resolution = 0.5,
        algorithm = 1
)

Idents(sub_5) <- "test_5"



p2 <- DimPlot(sub_5, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
pdf(file="cluster_5.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p2
dev.off()





DefaultAssay(sub_5) <- "RNA"
sub_5 <- NormalizeData(sub_5, normalization.method = "LogNormalize", scale.factor = 10000)

sub_5 <- FindVariableFeatures(sub_5, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(sub_5)
sub_5 <- ScaleData(sub_5, features = all.genes)


save(sub_5, file= "sub_RNA_5.RData")





