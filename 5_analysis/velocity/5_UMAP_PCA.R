#subset cluster

setwd("/analysis/velocity/subset_cluster")

#re-install seurat and sctransform
#remotes::install_version("Seurat", version = "4.0.5")
#remotes::install_version("sctransform", version = "0.3.3")


load("sub_RNA_5.RData")
library(Seurat)
library(SeuratObject)
library(purrr)
library(sctransform)

sub_5 <- RenameIdents(sub_5, 
                      
                      '0'='Suprabasal',
                      '1_0'='Secretory-1',
                      '1_1'='Secretory-1',
                      '1_2'='Secretory-1',
                      '1_3'='Secretory-1',
                      '1_4'='Secretory-1',
                      '1_5'='Goblet',
                      '1_6'='Secretory-1',
                      '1_7'='Goblet/Ionocyte',
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
                      '9_0'	="Ciliated+SARShi",
                      '9_1'	="Secretory-Ciliated+SARShi",
                      '10_0' ="Goblet+ISGhi",
                      '10_1' ="Secretory+ISGhi",
                      '10_2' ="Secretory+ISGhi",
                      '10_3' ="Secretory+ISGhi",
                      '10_4' ="Goblet+ISGhi",
                      '10_5' ="Secretory+ISGhi",
                      '10_6' ="Secretory+ISGhi",
                      '10_7' ="Secretory+ISGhi",
                      '11'	="Deuterosomal",
                      '12_0' = "Ionocyte",
                      '12_1' = "Brush/Tuft",
                      '13' = "Secretory-Ciliated"
)

sub_5@meta.data$cell_type <- as.character(Idents(sub_5))


##Subset data for only secretory, ciliated, secretory-ciliated

## S3 method for class 'Seurat'


sub_6 <- subset(
  sub_5,
  idents = c("Secretory-1", "Secretory-2", "Secretory+ISGhi", "Ciliated-1", "Ciliated-2", "Ciliated-3", "Ciliated+SARShi", "Secretory-Ciliated+SARShi", "Goblet+ISGhi",
             "Goblet", "Secretory-Ciliated"),
)


save(sub_6, file="sub_6.RData")

sub_6@meta.data$type_Sample <- paste0(sub_6@meta.data$type,"_", sub_6@meta.data$Sample)

list <- levels(as.factor(sub_6@meta.data$type_Sample))

#Subset merged data based on type and Sample

Long_read_UK_72hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_adult_Sample1_Sample1")    
Long_read_UK_72hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_adult_Sample2_Sample2")   
Long_read_UK_72hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_adult_Sample3_Sample3")  
Long_read_UK_72hpi_child_1 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_child_Sample4_Sample4")   
Long_read_UK_72hpi_child_2 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_child_Sample5_Sample5")   
Long_read_UK_72hpi_child_3 <- subset(sub_6, subset = type_Sample == "Long_read_UK_72hpi_child_Sample6_Sample6")  

Long_read_uninfected_adult_1 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_adult_Sample1_Sample1")
Long_read_uninfected_adult_2 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_adult_Sample2_Sample2")
Long_read_uninfected_adult_3 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_adult_Sample3_Sample3")
Long_read_uninfected_child_1 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_child_Sample4_Sample4")
Long_read_uninfected_child_2 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_child_Sample5_Sample5")
Long_read_uninfected_child_3 <- subset(sub_6, subset = type_Sample == "Long_read_uninfected_child_Sample6_Sample6")

Long_read_VIC01_48hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_adult_Sample1_Sample1")    
Long_read_VIC01_48hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_adult_Sample2_Sample2")   
Long_read_VIC01_48hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_adult_Sample3_Sample3")  
Long_read_VIC01_48hpi_child_1 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_child_Sample4_Sample4")   
Long_read_VIC01_48hpi_child_2 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_child_Sample5_Sample5")   
Long_read_VIC01_48hpi_child_3 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_48hpi_child_Sample6_Sample6")  

Long_read_VIC01_72hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_adult_Sample1_Sample1")    
Long_read_VIC01_72hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_adult_Sample2_Sample2")   
Long_read_VIC01_72hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_adult_Sample3_Sample3")  
Long_read_VIC01_72hpi_child_1 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_child_Sample4_Sample4")   
Long_read_VIC01_72hpi_child_2 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_child_Sample5_Sample5")   
Long_read_VIC01_72hpi_child_3 <- subset(sub_6, subset = type_Sample == "Long_read_VIC01_72hpi_child_Sample6_Sample6")  


Short_read_UK_72hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_adult_Sample1_Sample1")    
Short_read_UK_72hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_adult_Sample2_Sample2")   
Short_read_UK_72hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_adult_Sample3_Sample3")  
Short_read_UK_72hpi_child_1 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_child_Sample4_Sample4")   
Short_read_UK_72hpi_child_2 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_child_Sample5_Sample5")   
Short_read_UK_72hpi_child_3 <- subset(sub_6, subset = type_Sample == "Short_read_UK_72hpi_child_Sample6_Sample6")  

Short_read_uninfected_adult_1 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_adult_Sample1_Sample1")
Short_read_uninfected_adult_2 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_adult_Sample2_Sample2")
Short_read_uninfected_adult_3 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_adult_Sample3_Sample3")
Short_read_uninfected_child_1 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_child_Sample4_Sample4")
Short_read_uninfected_child_2 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_child_Sample5_Sample5")
Short_read_uninfected_child_3 <- subset(sub_6, subset = type_Sample == "Short_read_uninfected_child_Sample6_Sample6")

Short_read_VIC01_48hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_adult_Sample1_Sample1")    
Short_read_VIC01_48hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_adult_Sample2_Sample2")   
Short_read_VIC01_48hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_adult_Sample3_Sample3")  
Short_read_VIC01_48hpi_child_1 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_child_Sample4_Sample4")   
Short_read_VIC01_48hpi_child_2 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_child_Sample5_Sample5")   
Short_read_VIC01_48hpi_child_3 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_48hpi_child_Sample6_Sample6")  

Short_read_VIC01_72hpi_adult_1 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_adult_Sample1_Sample1")    
Short_read_VIC01_72hpi_adult_2 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_adult_Sample2_Sample2")   
Short_read_VIC01_72hpi_adult_3 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_adult_Sample3_Sample3")  
Short_read_VIC01_72hpi_child_1 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_child_Sample4_Sample4")   
Short_read_VIC01_72hpi_child_2 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_child_Sample5_Sample5")   
Short_read_VIC01_72hpi_child_3 <- subset(sub_6, subset = type_Sample == "Short_read_VIC01_72hpi_child_Sample6_Sample6")  


#Run sctransform on the list

list <- c(Long_read_UK_72hpi_adult_1, Long_read_UK_72hpi_adult_2, Long_read_UK_72hpi_adult_3,
          Long_read_UK_72hpi_child_1, Long_read_UK_72hpi_child_2, Long_read_UK_72hpi_child_3,
          Long_read_uninfected_adult_1, Long_read_uninfected_adult_2, Long_read_uninfected_adult_3,
          Long_read_uninfected_child_1, Long_read_uninfected_child_2, Long_read_uninfected_child_3,
          Long_read_VIC01_48hpi_adult_1, Long_read_VIC01_48hpi_adult_2, Long_read_VIC01_48hpi_adult_3,
          Long_read_VIC01_48hpi_child_1, Long_read_VIC01_48hpi_child_2, Long_read_VIC01_48hpi_child_3,
          Long_read_VIC01_72hpi_adult_1, Long_read_VIC01_72hpi_adult_2, Long_read_VIC01_72hpi_adult_3,
          Long_read_VIC01_72hpi_child_1, Long_read_VIC01_72hpi_child_2, Long_read_VIC01_72hpi_child_3,
          Short_read_UK_72hpi_adult_1, Short_read_UK_72hpi_adult_2, Short_read_UK_72hpi_adult_3,
          Short_read_UK_72hpi_child_1, Short_read_UK_72hpi_child_2, Short_read_UK_72hpi_child_3,
          Short_read_uninfected_adult_1, Short_read_uninfected_adult_2, Short_read_uninfected_adult_3,
          Short_read_uninfected_child_1, Short_read_uninfected_child_2, Short_read_uninfected_child_3,
          Short_read_VIC01_48hpi_adult_1, Short_read_VIC01_48hpi_adult_2, Short_read_VIC01_48hpi_adult_3,
          Short_read_VIC01_48hpi_child_1, Short_read_VIC01_48hpi_child_2, Short_read_VIC01_48hpi_child_3,
          Short_read_VIC01_72hpi_adult_1, Short_read_VIC01_72hpi_adult_2, Short_read_VIC01_72hpi_adult_3,
          Short_read_VIC01_72hpi_child_1, Short_read_VIC01_72hpi_child_2, Short_read_VIC01_72hpi_child_3)
          


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


fun1 <- function(x) {
  
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
  x <- RunPCA(x, features = VariableFeatures(x), ndims.print = 6:10, nfeatures.print = 10,approx=FALSE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  #Alternate Workflow regression
  x$CC.Difference <- x$S.Score - x$G2M.Score
  # run sctransform (with glmGamPoi)
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "CC.Difference", verbose = FALSE)
  x <- RunPCA(x, features = c(s.genes, g2m.genes),approx=FALSE)
  DimPlot(x)
  x

  
}

fin <- lapply(list, fun1)

save(fin, file = "sct_objects.RData")



#Rename each object

nm <- c("Long_read_UK_72hpi_adult_1", "Long_read_UK_72hpi_adult_2", "Long_read_UK_72hpi_adult_3",
          "Long_read_UK_72hpi_child_1", "Long_read_UK_72hpi_child_2", "Long_read_UK_72hpi_child_3",
          "Long_read_uninfected_adult_1", "Long_read_uninfected_adult_2", "Long_read_uninfected_adult_3",
          "Long_read_uninfected_child_1", "Long_read_uninfected_child_2", "Long_read_uninfected_child_3",
          "Long_read_VIC01_48hpi_adult_1", "Long_read_VIC01_48hpi_adult_2", "Long_read_VIC01_48hpi_adult_3",
          "Long_read_VIC01_48hpi_child_1", "Long_read_VIC01_48hpi_child_2", "Long_read_VIC01_48hpi_child_3",
         "Long_read_VIC01_72hpi_adult_1", "Long_read_VIC01_72hpi_adult_2", "Long_read_VIC01_72hpi_adult_3",
          "Long_read_VIC01_72hpi_child_1", "Long_read_VIC01_72hpi_child_2", "Long_read_VIC01_72hpi_child_3",
          "Short_read_UK_72hpi_adult_1", "Short_read_UK_72hpi_adult_2", "Short_read_UK_72hpi_adult_3",
          "Short_read_UK_72hpi_child_1", "Short_read_UK_72hpi_child_2", "Short_read_UK_72hpi_child_3",
          "Short_read_uninfected_adult_1", "Short_read_uninfected_adult_2", "Short_read_uninfected_adult_3",
          "Short_read_uninfected_child_1", "Short_read_uninfected_child_2", "Short_read_uninfected_child_3",
          "Short_read_VIC01_48hpi_adult_1", "Short_read_VIC01_48hpi_adult_2", "Short_read_VIC01_48hpi_adult_3",
          "Short_read_VIC01_48hpi_child_1", "Short_read_VIC01_48hpi_child_2", "Short_read_VIC01_48hpi_child_3",
         "Short_read_VIC01_72hpi_adult_1", "Short_read_VIC01_72hpi_adult_2", "Short_read_VIC01_72hpi_adult_3",
          "Short_read_VIC01_72hpi_child_1", "Short_read_VIC01_72hpi_child_2", "Short_read_VIC01_72hpi_child_3")



names(fin) <- nm


# Isolate the datasets in nested list

Long_read_UK_72hpi_adult_1 <-  fin[[1]]
Long_read_UK_72hpi_adult_2 <-  fin[[2]]
Long_read_UK_72hpi_adult_3 <- fin[[3]]
Long_read_UK_72hpi_child_1 <- fin[[4]] 
Long_read_UK_72hpi_child_2 <-  fin[[5]]
Long_read_UK_72hpi_child_3 <- fin[[6]]
Long_read_uninfected_adult_1 <-  fin[[7]]
Long_read_uninfected_adult_2 <-  fin[[8]]
Long_read_uninfected_adult_3 <- fin[[9]]
Long_read_uninfected_child_1 <-  fin[[10]]
Long_read_uninfected_child_2 <-  fin[[11]]
Long_read_uninfected_child_3 <- fin[[12]]
Long_read_VIC01_48hpi_adult_1 <-  fin[[13]]
Long_read_VIC01_48hpi_adult_2 <-  fin[[14]]
Long_read_VIC01_48hpi_adult_3 <- fin[[15]]
Long_read_VIC01_48hpi_child_1 <-  fin[[16]]
Long_read_VIC01_48hpi_child_2 <-  fin[[17]]
Long_read_VIC01_48hpi_child_3 <- fin[[18]]
Long_read_VIC01_72hpi_adult_1 <-  fin[[19]]
Long_read_VIC01_72hpi_adult_2 <-  fin[[20]]
Long_read_VIC01_72hpi_adult_3 <- fin[[21]]
Long_read_VIC01_72hpi_child_1 <- fin[[22]] 
Long_read_VIC01_72hpi_child_2 <-  fin[[23]]
Long_read_VIC01_72hpi_child_3 <- fin[[24]]
Short_read_UK_72hpi_adult_1 <-  fin[[25]]
Short_read_UK_72hpi_adult_2 <-  fin[[26]]
Short_read_UK_72hpi_adult_3 <- fin[[27]]
Short_read_UK_72hpi_child_1 <-  fin[[28]]
Short_read_UK_72hpi_child_2 <-  fin[[29]]
Short_read_UK_72hpi_child_3 <- fin[[30]]
Short_read_uninfected_adult_1 <-  fin[[31]]
Short_read_uninfected_adult_2 <-  fin[[32]]
Short_read_uninfected_adult_3 <- fin[[33]]
Short_read_uninfected_child_1 <-  fin[[34]]
Short_read_uninfected_child_2 <-  fin[[35]]
Short_read_uninfected_child_3 <- fin[[36]]
Short_read_VIC01_48hpi_adult_1 <-  fin[[37]]
Short_read_VIC01_48hpi_adult_2 <-  fin[[38]]
Short_read_VIC01_48hpi_adult_3 <- fin[[39]]
Short_read_VIC01_48hpi_child_1 <-  fin[[40]]
Short_read_VIC01_48hpi_child_2 <-  fin[[41]]
Short_read_VIC01_48hpi_child_3 <- fin[[42]]
Short_read_VIC01_72hpi_adult_1 <-  fin[[43]]
Short_read_VIC01_72hpi_adult_2 <-  fin[[44]]
Short_read_VIC01_72hpi_adult_3 <- fin[[45]]
Short_read_VIC01_72hpi_child_1 <-  fin[[46]]
Short_read_VIC01_72hpi_child_2 <-  fin[[47]]
Short_read_VIC01_72hpi_child_3 <- fin[[48]]




##merge the data


fin_merge <-  merge(Long_read_UK_72hpi_adult_1, c(Long_read_UK_72hpi_adult_2,Long_read_UK_72hpi_adult_3,
                                                                                  Long_read_uninfected_adult_1,Long_read_uninfected_adult_2,Long_read_uninfected_adult_3,
                                                                                  Long_read_VIC01_48hpi_adult_1,Long_read_VIC01_48hpi_adult_2,Long_read_VIC01_48hpi_adult_3,
                                                                                  Long_read_VIC01_72hpi_adult_1, Long_read_VIC01_72hpi_adult_2, Long_read_VIC01_72hpi_adult_3,
                                                                                  Short_read_UK_72hpi_adult_1,  Short_read_UK_72hpi_adult_2,  Short_read_UK_72hpi_adult_3,
                                                                                  Short_read_uninfected_adult_1,  Short_read_uninfected_adult_2,  Short_read_uninfected_adult_3,
                                                                                  Short_read_VIC01_48hpi_adult_1,  Short_read_VIC01_48hpi_adult_2,  Short_read_VIC01_48hpi_adult_3,
                                                                                  Short_read_VIC01_72hpi_adult_1,  Short_read_VIC01_72hpi_adult_2,  Short_read_VIC01_72hpi_adult_3,
                                                                                  Long_read_UK_72hpi_child_1, Long_read_UK_72hpi_child_2,Long_read_UK_72hpi_child_3,
                                                                                  Long_read_uninfected_child_1,Long_read_uninfected_child_2,Long_read_uninfected_child_3,
                                                                                  Long_read_VIC01_48hpi_child_1,Long_read_VIC01_48hpi_child_2,Long_read_VIC01_48hpi_child_3,
                                                                                  Long_read_VIC01_72hpi_child_1, Long_read_VIC01_72hpi_child_2, Long_read_VIC01_72hpi_child_3,
                                                                                  Short_read_UK_72hpi_child_1,  Short_read_UK_72hpi_child_2,  Short_read_UK_72hpi_child_3,
                                                                                  Short_read_uninfected_child_1,  Short_read_uninfected_child_2,  Short_read_uninfected_child_3,
                                                                                  Short_read_VIC01_48hpi_child_1,  Short_read_VIC01_48hpi_child_2,  Short_read_VIC01_48hpi_child_3,
                                                                                  Short_read_VIC01_72hpi_child_1,  Short_read_VIC01_72hpi_child_2,  Short_read_VIC01_72hpi_child_3))


#fin_merge <- merge(fin)

##https://github.com/satijalab/seurat/issues/2814

VariableFeatures(fin_merge[["SCT"]]) <- rownames(fin_merge[["SCT"]]@scale.data)


save(fin_merge, file = "fin_merge.RData")




fin_merge<- RunPCA(fin_merge,verbose = FALSE, approx=FALSE)
fin_merge<- RunUMAP(fin_merge, reduction = "pca", dims = 1:10, verbose = FALSE)
fin_merge<- FindNeighbors(fin_merge, reduction = "pca", dims = 1:10, verbose = FALSE)
fin_merge<- FindClusters(fin_merge, resolution = 0.3, verbose = FALSE)

ElbowPlot(fin_merge)

# Visualization
p1 <- DimPlot(fin_merge, reduction = "umap", group.by = "type", raster=TRUE)
pdf(file="type.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p1
dev.off()

# Visualization
p2 <- DimPlot(fin_merge, reduction = "umap", label = TRUE, repel = FALSE, raster=FALSE)
pdf(file="cluster.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p2
dev.off()

p9 <- DimPlot(fin_merge, reduction = "umap", group.by = "Infection", raster=TRUE)
pdf(file="Infection.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p9
dev.off()

p11 <- FeaturePlot(fin_merge, reduction = "umap", features = "percent_mito", raster=TRUE)
pdf(file="mito.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p11
dev.off()

p12 <- DimPlot(fin_merge, reduction = "umap", group.by = "Phase")
pdf(file="ccphase.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p12
dev.off()

p13 <- DimPlot(fin_merge, reduction = "umap", group.by = "Infection_tier", raster=TRUE, cols = c("grey", "#f0f9e8", "#7bccc4","#43a2ca","#0868ac")) +
  theme(axis.text= element_text(size=22, color="black", face="bold"), axis.title=element_text(size=22, color="black", face="bold"),
        legend.text = element_text(size=10, color="black", face="bold"), text=element_text(size=22, color="black", face="bold"))

p13



save(fin_merge, file = "fin_merge_post_UMAP.RData")

##Subcluster 8
DefaultAssay(fin_merge) <- "SCT"


c2 <- FindSubCluster(
  fin_merge,
  "8",
  graph.name = "SCT_snn",
  subcluster.name = "sp_8",
  resolution = 0.2,
  algorithm = 1
)
Idents(c2) <- c2@meta.data$sp_8

DimPlot(c2, reduction = "umap", label=TRUE)



levels(c2) <- sort(levels(c2), decreasing=FALSE)

save(c2, file = "recluster_sub_dims10_res_0.3_c2.RData")



##Subcluster 8
DefaultAssay(fin_merge) <- "SCT"


c3 <- FindSubCluster(
  fin_merge,
  "8",
  graph.name = "SCT_snn",
  subcluster.name = "sp_8_2",
  resolution = 0.3,
  algorithm = 1
)
Idents(c3) <- c3@meta.data$sp_8_2

DimPlot(c3, reduction = "umap", label=TRUE)



levels(c3) <- sort(levels(c3), decreasing=FALSE)

save(c3, file = "recluster_sub_dims10_res_0.3_c3.RData")


DefaultAssay(c3) <- "RNA"
c3 <- NormalizeData(c3, normalization.method = "LogNormalize", scale.factor = 10000)

c3 <- FindVariableFeatures(c3, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(c3)
c3 <- ScaleData(c3, features = all.genes)

save(c3, file= "c3_RNA.RData")



DefaultAssay(c3) <- "RNA"

p1 <- DotPlot(c3, features=c("KRT5", "TP63", "MKI67", "KRT13", "KRT4","FOXJ1", "SNTN", "CAPS", "SCGB1A1", "MUC5AC","MUC5B", "DEUP1", "LRMP", "ASCL3", "CFTR")) &coord_flip()
p1+theme(axis.title = element_text(face=,"bold", size=17), axis.text=element_text(face="bold", size=17), text = element_text(face="bold", size=15),
         legend.text=element_text(face="bold", size=15), legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))

p2 <- FeaturePlot(c3, features=c("FOXJ1", "SNTN", "CAPS", "SCGB1A1", "MUC5AC","MUC5B"))
p2



##Analysis


#Findmarkers 5 vs 6


f1 <- FindMarkers(c3, ident.1 = "5", ident.2= "6", min.pct=0.25, min.diff.pct=0.25)
#Warning message:
 # In FindMarkers.default(object = data.use, slot = data.slot, counts = counts,  :
  #No features pass min.diff.pct threshold; returning empty data.frame


f1 <- FindMarkers(c3, ident.1 = "5", ident.2= "6", min.pct=0.25)
f1

genes <- rownames(f1)


DefaultAssay(c3) <- "RNA"

c3.markers <- FindAllMarkers(c3, only.pos = TRUE)
save(c3.markers, file="c3.markers.RData")



##Set clusters


DefaultAssay(c3) <- "RNA"
c3 <- RenameIdents(c3, 
                      
                      '0'='Secretory',
                      '1'= 'Ciliated',
                      '2'	="Secretory",
                      '3'	="Ciliated",
                      '4'	="Secretory",
                      '5'	="Secretory-Ciliated",
                      '6'	="Secretory-Ciliated",
                      '7'	="Secretory+ISGhi",
                      '8_0'	="Ciliated+SARShi",
                      '8_1'	="Ciliated+SARShi",
                      '8_2'	="Secretory-Cilited+SARShi",
                      '8_3'	="Ciliated+SARShi",
                      '8_4'	="Ciliated+SARShi"
                      
)

c3@meta.data$cell_type_recluster <- as.character(Idents(c3))

levels(c3)

save(c3, file="c3_RNA_renamed.RData")


VlnPlot(c3, features=c("FOXJ1", "SNTN", "CAPS", "SCGB1A1", "MUC5AC","MUC5B"), pt.size = 0)
DotPlot(c3, features=c("FOXJ1", "SNTN", "CAPS", "SCGB1A1", "MUC5AC","MUC5B"))&coord_flip()
FeaturePlot(c3,features=c("FOXJ1", "SNTN", "CAPS")) +p2

AverageExpression(c3,features=c("FOXJ1", "SNTN", "CAPS", "SCGB1A1", "MUC5AC","MUC5B" ))

DefaultAssay(c3) <- "SCT"
# Visualization

p2 <- DimPlot(c3, reduction = "umap", label = TRUE, repel = FALSE, raster=FALSE)
pdf(file="dimplot_recluster_c3.pdf",  width = 11.69, # The width of the plot in inches
    height = 8.27)
p2
dev.off()
