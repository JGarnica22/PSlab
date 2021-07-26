library(SeuratWrappers)
library(monocle3)

# Use seurat wrappers to convert your Seurat object where you can label cells as you want (cell type, kmeams...) into an
# cds object required for monocle analysis

cds <- as.cell_data_set(scdata)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, verbose = T)
plot_pc_variance_explained(cds)
# Compute gene expression
cds <- estimate_size_factors(cds, round_exprs = F)
rowData(cds)$gene_short_name <- row.names(rowData(cds))

## Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
seur_clus <- plot_cells(cds, color_cells_by = "celltype",
                        label_cell_groups = F, cell_size = 1)
                        
## Step 4: Cluster cells on monocle where to base the trajectory analysis                
set.seed(4)
cds <- cluster_cells(cds, resolution = 2.5e-5)
mono_clusters <- plot_cells(cds, color_cells_by = "cluster",
                            label_cell_groups = F, cell_size = 1)

## Step 5: Peform trajectory analysis
cds <- learn_graph(cds)
cds <- order_cells(cds)
# Plot trajectoryt with the grouping you like
um <- plot_cells(cds, color_cells_by = "celltype",
                 label_cell_groups = F, cell_size = 1)
# Plot cells with pseudotime values
um_ps <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1)



diff_test_res <- differentialGeneTest(cds,
                    fullModelFormulaStr = "~sm.ns(Pseudotime)")

## Loop to generate all the figures in a pdf for a list of genes inclucded in "features"
pdf("figs/monocle3_analysis.pdf", width = 11, height = 11)
mono_clusters
seur_clus
um
um_ps
for (i in seq(1, 106, 5)){
  if (i+4 < length(features)){
  table_cds <- cds[rowData(cds)$gene_short_name %in%
                  features[i:c(i+4)],]
  } else {table_cds <- cds[rowData(cds)$gene_short_name %in%
                             features[c(length(features)-4):length(features)],]}
  fp <- plot_genes_in_pseudotime(table_cds, color_cells_by = "celltype")+
    ggplot2::ggtitle("Gene expression over pseudotime")
  print(fp)
}
for (i in seq(1, 106, 9)){
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
