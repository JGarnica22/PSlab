library(SeuratWrappers)
library(monocle3)
#conver Seurat object into single cell dataset:
cds <- as.cell_data_set(scdata)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, verbose = T)
plot_pc_variance_explained(cds)

# estimate gene expression data
cds <- estimate_size_factors(cds, round_exprs = F)
rowData(cds)$gene_short_name <- row.names(rowData(cds))

## Step 2(optional): Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "batch")
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

#Plot cells based on Seurat classification
plot_cells(cds, color_cells_by = "sampletype",
                        label_cell_groups = F, cell_size = 1) +
  ggtitle("Celltype based on 3 kmeans")

# Plot a featureplot to assess a correct grouping
plot_cells(cds, genes= "Il10")

# Cluster cells to generate the trajectory
# We need to try different resolution and seeding to obtain a similar
# clustering to Seurat.
r <- 6e-3
# for (i in 1:24){
  s <- 1
  set.seed(s)
  cds <- cluster_cells(cds, resolution = r)
  cds <- learn_graph(cds)
  mono_clusters <- plot_cells(cds, color_cells_by = "cluster",
                              label_cell_groups = F, cell_size = 1) +
    ggtitle(paste0("Monocle3 clustering. Seed: ", s, ". Res = ", r,". ",
                   clu))
  print(mono_clusters)
#}

# Decide where to start the trajectory
cds <- order_cells(cds)

# Plot results as UMAPs
um <- plot_cells(cds, color_cells_by = "sampletype",
                 label_cell_groups = F, cell_size = 1) +
  ggplot2::ggtitle(paste0("Only ", clu, " trajectory. ",
                          "Seed: ", s, ". Res = ", r))
um_ps <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1) +
  ggplot2::ggtitle(paste0("Only ", clu, " pseudotime. ",
                          "Seed: ", s, ". Res = ", r))


# Plot feature plots on UMAP and gene expression over pseudotime
# for all the genes in 'features' list.
pdf(paste0("figs/monocle3_analysis_", clu, ".pdf"),
    width = 11, height = 11)
um
um_ps
for (i in seq(1, 106, 5)){
  if (i+4 < length(features)){
    table_cds <- cds[rowData(cds)$gene_short_name %in%
                       features[i:c(i+4)],]
  } else {table_cds <- cds[rowData(cds)$gene_short_name %in%
                             features[c(length(features)-4):length(features)],]}
  fp <- plot_genes_in_pseudotime(table_cds, color_cells_by = "sampletype")+
    ggplot2::ggtitle("Gene expression over pseudotime")
  print(fp)
}
for (i in seq(1, length(features), 9)){
  if (i+8 < length(features)){
    f <- plot_cells(cds,
                    genes= features[i:c(i+8)]) +
      ggplot2::ggtitle("Gene expression features plots")
  } else {f <- plot_cells(cds,
                          genes= features[c(length(features)-8):length(features)]) +
    ggplot2::ggtitle("Gene expression features plots")}
  print(f)
}
dev.off()

# Generate table with Morans' I value computing the pseudotime
diff_5 <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
writexl::write_xlsx(diff_5, "out/diff_pseudotimel_5.xlsx")
