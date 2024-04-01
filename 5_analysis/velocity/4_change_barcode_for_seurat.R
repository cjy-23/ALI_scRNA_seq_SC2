setwd("/velocity/velocyto_filtered/")
library(Seurat)

load("/velocity/velocyto_filtered/merged_subset.RData")
sp <- split(merged_ob@meta.data, merged_ob@meta.data$Type)

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

#Bind all nested data frames
meta_new <- do.call(rbind,sp)

nm <- names(sp)
nm


#Set as new metadata

merged_ob@meta.data <- meta_new

save(merged_ob, file="merged_subset_fix_barcode.RData")
