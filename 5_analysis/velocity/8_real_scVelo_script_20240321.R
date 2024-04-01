#R 4.3.2, Seurat_5.0.2, SeuratObject_5.0.1

setwd("/velocity/real_scVelo_20240321/")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(reticulate)
library(hdf5r)
library(loomR)

#devtools::install_github(repo = "hhoeflin/hdf5r")
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

load("/velocity/scVelo_run_20240314/merged_ob4_10000_20240314.RData")

levels(Idents(merged_ob4))
merged_ob4 <- RenameIdents(merged_ob4, "Secretory-Cilited+SARShi"= "Secretory-Ciliated+SARShi")

levels(Idents(merged_ob4))

merged_ob4@meta.data$cell_type_recluster <- Idents(merged_ob4)

levels(merged_ob4@meta.data$cell_type_recluster)

save(merged_ob4, file="/velocity/real_scVelo/merged_ob4_10000_20240321_revised.RData")

# Convert to loom file
loomR-package::write.loom(merged_ob4, filename = "merged_ob4_10000_20240320.loom", assays = c("spliced", "unspliced", "ambiguous"))
pfile <- SeuratDisk::Convert(merged_ob4, to = "loom", filename = "merged_ob4_10000_20240320.loom", 
                 display.progress = FALSE,  assay = c("spliced", "unspliced", "ambiguous"))



loom <- h5read("path/to/your/file.loom", "/")


# assuming that you have some Seurat object called merged_ob4:

# save metadata table:
merged_ob4$barcode <- colnames(merged_ob4)
merged_ob4$UMAP_1 <- merged_ob4@reductions$umap@cell.embeddings[,1]
merged_ob4$UMAP_2 <- merged_ob4@reductions$umap@cell.embeddings[,2]
write.csv(merged_ob4@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix ?? will change to spliced if only the gene names are derived from this
library(Matrix)
counts_matrix <- GetAssayData(merged_ob4, assay='spliced', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(merged_ob4@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)


write.csv(x = t(as.matrix(merged_ob4@assays$spliced@counts)), file = 'E11_spliced.csv')
write.csv(x = t(as.matrix(merged_ob4@assays$unspliced@counts)), file = 'E11_unspliced.csv')


