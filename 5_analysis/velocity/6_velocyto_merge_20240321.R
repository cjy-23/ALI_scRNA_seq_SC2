setwd("/velocity/velocyto_umap_final/")

library(Seurat)
library(dplyr)
##WITH R4.1.2 SeuratObject_4.1.3 Seurat_4.3.0  
load("/velocity/subset_cluster/dim10_res_0.3/c3_RNA_renamed.RData") #clustering
load("/velocity/velocyto_filtered/merged_subset_fix_barcode.RData") #velocyto loom file --> seurat



#change barcodes to match veloctyo loom files
meta <- c3@meta.data
meta$barcode <- rownames(meta)
meta$barcode_af1 <- gsub("-.*", "", meta$barcode)


sp <- split(meta, meta$type)


sp[[1]]$barcode_af2 <- paste0(sp[[1]]$barcode_af1, "_", "1")
sp[[2]]$barcode_af2 <- paste0(sp[[2]]$barcode_af1, "_", "2")
sp[[3]]$barcode_af2 <- paste0(sp[[3]]$barcode_af1, "_", "3")
sp[[4]]$barcode_af2 <- paste0(sp[[4]]$barcode_af1, "_", "4")
sp[[5]]$barcode_af2 <- paste0(sp[[5]]$barcode_af1, "_", "5")
sp[[6]]$barcode_af2 <- paste0(sp[[6]]$barcode_af1, "_", "6")
sp[[7]]$barcode_af2 <- paste0(sp[[7]]$barcode_af1, "_", "7")
sp[[8]]$barcode_af2 <- paste0(sp[[8]]$barcode_af1, "_", "8")
sp[[9]]$barcode_af2 <- paste0(sp[[9]]$barcode_af1, "_", "9")
sp[[10]]$barcode_af2 <- paste0(sp[[10]]$barcode_af1, "_", "10")
sp[[11]]$barcode_af2 <- paste0(sp[[11]]$barcode_af1, "_", "11")
sp[[12]]$barcode_af2 <- paste0(sp[[12]]$barcode_af1, "_", "12")
sp[[13]]$barcode_af2 <- paste0(sp[[13]]$barcode_af1, "_", "13")
sp[[14]]$barcode_af2 <- paste0(sp[[14]]$barcode_af1, "_", "14")
sp[[15]]$barcode_af2 <- paste0(sp[[15]]$barcode_af1, "_", "15")
sp[[16]]$barcode_af2 <- paste0(sp[[16]]$barcode_af1, "_", "16")

names(sp)

##Reorder nested list based on order in c3

sp2 <- list(sp$Long_read_UK_72hpi_adult, sp$Long_read_uninfected_adult, sp$Long_read_VIC01_48hpi_adult, sp$Long_read_VIC01_72hpi_adult,
            sp$Short_read_UK_72hpi_adult, sp$Short_read_uninfected_adult, sp$Short_read_VIC01_48hpi_adult, sp$Short_read_VIC01_72hpi_adult,
            sp$Long_read_UK_72hpi_child, sp$Long_read_uninfected_child, sp$Long_read_VIC01_48hpi_child, sp$Long_read_VIC01_72hpi_child,
            sp$Short_read_UK_72hpi_child, sp$Short_read_uninfected_child, sp$Short_read_VIC01_48hpi_child, sp$Short_read_VIC01_72hpi_child)


nm2 <- c("Long_read_UK_72hpi_adult", "Long_read_uninfected_adult", "Long_read_VIC01_48hpi_adult", "Long_read_VIC01_72hpi_adult",
            "Short_read_UK_72hpi_adult", "Short_read_uninfected_adult", "Short_read_VIC01_48hpi_adult", "Short_read_VIC01_72hpi_adult",
            "Long_read_UK_72hpi_child", "Long_read_uninfected_child", "Long_read_VIC01_48hpi_child", "Long_read_VIC01_72hpi_child",
            "Short_read_UK_72hpi_child", "Short_read_uninfected_child", "Short_read_VIC01_48hpi_child", "Short_read_VIC01_72hpi_child")

names(sp2) <- nm2

#Bind all nested data frames
meta_new <- do.call(rbind,sp2)


#Set as new metadata

c3@meta.data <- meta_new

##Check if we have same type cell barcode suffix (correct)

c3_meta <- c3@meta.data
merged_ob_meta <- merged_ob@meta.data

c3_meta_sp <- split(c3_meta, c3_meta$type)
merged_ob_sp <- split(merged_ob_meta, merged_ob_meta$Type)


##Just double check whether we have same barcodes

c3_names <- c3@meta.data$barcode_af2
merged_ob_names <- merged_ob@meta.data$barcode_af2

c3_names_df <- as.data.frame(c3_names)

setdiff(c3_names, merged_ob_names) ##fine



#Rename the cells

c3_r <- RenameCells(c3, new.names = c3_names)

merged_ob2 <- RenameCells(merged_ob, new.names = merged_ob_names)


identical(c3@meta.data$cell_type_recluster, c3_r@meta.data$cell_type_recluster) ##TRUE

c3_r@meta.data$cell_type_recluster <- as.character(c3_r@meta.data$cell_type_recluster)
DimPlot(c3_r, group.by="cell_type_recluster")


identical(rownames(c3_r@meta.data), c3_names)
save(c3_r, file="c3_RNA_renamed_barcode_20240314.RData")


##Make sure the cell barcode order is the same between two seurat objects
barcodes1_2 <- colnames(c3)
barcodes2_2 <- colnames(merged_ob)
are_equal_order <- identical(barcodes1_2, barcodes2_2) ##false!

barcodes1_3 <- colnames(c3_r)
barcodes2_3 <- colnames(merged_ob2)

are_equal_order <- identical(barcodes1_3, barcodes2_3) ##false!


##======================================
#Add PCA & umap embeddings to veloctyo seurat object

merged_ob2@reductions[["umap"]] <- c3_r@reductions[["umap"]]
merged_ob2@reductions[["pca"]] <- c3_r@reductions[["pca"]]


DimPlot(c3_r, group.by="cell_type_recluster") ##Fine


ggg <- merged_ob2@meta.data

matching_indices <- match(ggg$barcode_af2, c3_r@meta.data$barcode_af2)
ggg$barcode_af2[1]
c3_r@meta.data$barcode_af2[23]
matched_cell_types <- c3_r@meta.data$cell_type_recluster[matching_indices]

# Assign matched cell types where barcodes match
ggg$cell_type_recluster[!is.na(matching_indices)] <- matched_cell_types[!is.na(matched_cell_types)]


merged_ob2@meta.data <- ggg

DimPlot(merged_ob2, group.by="cell_type_recluster")

Idents(merged_ob2) <- merged_ob2@meta.data$cell_type_recluster
DimPlot(merged_ob2)

save(merged_ob2, file="final_merged_velocity_umap_20240314.RData")
