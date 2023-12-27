
library(monocle3)
library(SeuratWrappers)


##Install proper versions of packages:

#remotes::install_version("SeuratObject", version = "4.1.2")
#remotes::install_version("Seurat", version = "4.2.0")

##Load data
small_monocle <- readRDS("small_monocle.Rds")


plot_cells(small_monocle,
           color_cells_by = "celltype",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
small_monocle <- order_cells(small_monocle)

pdf("monocle.pdf",7,6)
plot_cells(small_monocle,
           color_cells_by = "pseudotime",
           group_cells_by = "celltype",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "black")
dev.off()


pdf("leaf_monocle.pdf",7,6)
plot_cells(small_monocle,
           color_cells_by = "pseudotime",
           group_cells_by = "celltype",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "black")
dev.off()

