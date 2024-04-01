load("sub_RNA_5.RData")


library(Seurat)
library(tidyverse)

dim(sub_5)

small = sub_5[, sample(123128, 10000, replace = FALSE)]

small$celltype <- recode(small@active.ident,

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
                      '9_0'	="Ciliated+SC2",
                      '9_1'	="Secretory-Ciliated+SC2",
                      '10_0' ="Goblet+IFN-stim",
                      '10_1' ="Secretory-3",
                      '10_2' ="Secretory-3",
                      '10_3' ="Secretory-3",
                      '10_4' ="Goblet+IFN-stim",
                      '10_5' ="Secretory-3",
                      '10_6' ="Secretory-3",
                      '10_7' ="Secretory-3",
                      '11'	="Deuterosomal",
                      '12_0' = "Ionocyte",
                      '12_1' = "Brush/Tuft",
                      '13' = "Secretory-Ciliated"
)


library(monocle3)
library(SeuratWrappers)


cds <- as.cell_data_set(small)

plot_cells(cds, color_cells_by = "celltype", show_trajectory_graph = FALSE)

cds <- cluster_cells(cds, resolution=1e-3)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
```
```{r}
cds <- order_cells(cds)

pdf("../plots/monocle.pdf",7,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "celltype",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")
dev.off()


##Load data
#small_monocle <- readRDS("small_monocle.Rds")


plot_cells(small_monocle,
           color_cells_by = "celltype",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
small_monocle <- order_cells(small_monocle)



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

