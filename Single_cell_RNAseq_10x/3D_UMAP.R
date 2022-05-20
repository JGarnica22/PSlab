#### 
## Date: 20 - 05 - 2022
###

###
## Description: 3D Plot of Seurat clusters using monocle3. It assumes the Seurat object has gone through
## the previous processes (qc, transformations and clustering). Reference used in this example is available
## in the data section of the repository
###

library(Seurat)
library(monocle3)
library(scales)

#Read seurat object 
seurat.object <- readRDS("data/vision_reference.Rds")
DefaultAssay(seurat.object) <- "RNA"

#Transform seurat object to cell datsa set for monocle3
cds <- SeuratWrappers::as.cell_data_set(seurat.object)
#Preprocess
cds <- preprocess_cds(cds, num_dim = 20, verbose = T)
cds <- estimate_size_factors(cds, round_exprs = F)
rowData(cds)$gene_short_name <- row.names(rowData(cds))

#Align to minimize batch caused for using different samples
cds <- align_cds(cds, num_dim = 10, alignment_group = "sampletype")

#Dimensional reduction and cluster prediction
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d, resolution = 0.0002)

#Ggplot-like palette. Using as many colors as unique clusters. In this object, clusters are stored in
#vision meta data column
pal <- hue_pal()(9)
cds_3d <- learn_graph(cds_3d, use_partition = F)
cds_3d <- order_cells(cds_3d)
pseudo_3d <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")
clusters_3d <- plot_cells_3d(cds_3d, color_cells_by = "vision", show_trajectory_graph = F, 
                          color_palette = pal)
