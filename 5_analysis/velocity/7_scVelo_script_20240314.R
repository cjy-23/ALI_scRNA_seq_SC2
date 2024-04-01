#Sys.setenv(PKG_CXXFLAGS="-DARMA_64BIT_WORD")

#R 4.3.2, Seurat_5.0.2, SeuratObject_5.0.1


setwd("/velocity/scVelo_run_20240314/")

library(velocyto.R)
library(SeuratWrappers)
library(Seurat)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
#install.packages("RcppArmadillo", repos=NULL)
#install.packages("Rcpp")
#install.packages("R.utils")
#remotes::install_github('satijalab/seurat-wrappers')
#install_github("velocyto-team/velocyto.R")
#BiocManager::install("pcaMethods")
#install.packages("RcppArmadillo")

load("/velocity/velocyto_umap_final/final_merged_velocity_umap_20240314.RData")


#DefaultAssay(merged_ob2) <- "SCT"
merged_ob4 = merged_ob2[, sample(71818, 10000, replace = FALSE)]

save(merged_ob4, file="merged_ob4_10000_20240314.RData")


merged_ob5 = RunVelocity(object = merged_ob4, deltaT = 1, kCells = 10, fit.quantile = 0.02 )


save(merged_ob5, file="merged_ob5_velo.RData")

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = merged_ob5)))
names(x = ident.colors) <- levels(x = merged_ob5)
cell.colors <- ident.colors[Idents(object = merged_ob5)]
names(x = cell.colors) <- colnames(x = merged_ob5)

vel_input <- Tool(object = merged_ob5, slot = "RunVelocity")

pdf("velo.pdf")
show.velocity.on.embedding.cor(emb = Embeddings(object = merged_ob5, reduction = "umap"),
                               vel = vel_input,
                               n = 50, scale = "log", cell.colors = ac(x = cell.colors, alpha = 0.5),
                               cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40,
                               arrow.lwd = 1, cell.border.alpha = 0.1) 
dev.off()
