#This script is to analyze Single-cell RNAseq data from aggr (Cell Ranger output)
#It was created for R 3.6.3 version (2020-05-22)
#Copyright (C) 2020  Patricia Solé Sánchez
#################################################################################################

# Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. 
# Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell 
# transcriptomic measurements, and to integrate diverse types of single-cell data.

# https://satijalab.org/seurat/

cran.packages <- c("Seurat", "dplyr", "patchwork", "devtools", "cowplot", "magrittr")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

bioc.packages <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'limma', 'S4Vectors', 'SingleCellExperiment',
                   'SummarizedExperiment', 'batchelor', 'Matrix.utils')
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

library(devtools)

dev.pack <- c("satijalab/seurat-data")
for (i in dev.pack) {
  if(!require(i, character.only = TRUE)) {
    devtools::install_github(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

library(Seurat)
library(SeuratData)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(magrittr)

# Set your working directory (the project you are working in):
setwd("/Users/patri/Documents/LAB/TESIS DOCTORAL 2015-2020/TR1 PROJECT/2020_05_10xGenomics_BDC_INS13-21_SANTAMARIA_09/Seurat")

# Specify parameters to be used along the script:
#Indicate species working with (typicall mouse or human)
species <- "mouse"
#Indicate samples you sequenced separatedly (for example, you sorted your treated
#and control populations separatedly):
#Indicate samples in the order they are barcoded!!
#It's better to indicate sample type as TREATMENT_CONDITION
sampletype <- c("BDC_CTL", "BDC_TET", "INS_CTL", "INS_TET")
#Usually, those samples can be identified by the cell name. Each cell has a -number
#at the end that identifies the sample.

# Read the data:
# The Read10X function reads in the output of the cellranger pipeline from 10X 
# (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
# returning a unique molecular identified (UMI) count matrix. The values in this 
# matrix represent the number of molecules for each feature (i.e. gene; row) that 
# are detected in each cell (column).
# Use filtered data

# Indicate in the Read10X function, the directory of your input files (3 files must be
# found in the data.dir: matrix.mt.gz, barcodes.tsv.gz and features.tsv.gz)
data <- Read10X(data.dir = paste0(getwd(), "/data/aggr/filtered_feature_bc_matrix/"))

# We next use the count matrix to create a Seurat object. The object serves as a 
# container that contains both data (like the count matrix) and analysis (like PCA, 
# or clustering results) for a single-cell dataset.
# min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. 
# min.features: Include cells where at least this many features are detected.
scdata <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)

# Check the expression (counts) for some genes in the first 100 cells:
# Change genes for those you want to check:
data[c("Sell", "Cd3e", "Il10"), 1:100]

# The . values in the matrix represent 0s (no molecules detected). Since most values in 
# an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. 
# This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

# Standard pre-processing workflow for scRNA-seq data in Seurat:
# 1) Selection and filtration of cells based on QC metrics
# 2) Data normalization and scaling
# 3) Detection of highly variable features

# QC and selecting cells for further analysis
############################################################
#Criteria:
# * The number of unique genes detected in each cell (Low-quality cells or empty droplets will 
# often have very few genes and cell doublets or multiplets may exhibit an aberrantly high gene count
# * The percentage of reads that map to the mitochondrial genome (Low-quality/dying cells often
# exhibit extensive mitochondrial contamination

# Calculate mitochondrial QC metrics with the PercentageFeatureSet function:
if (species == "mouse") {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^mt-")
} else {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
}

# QC metrics can be found in the Seurat object metadata:
head(scdata@meta.data)

# Visualize QC metrics as a violin plot:
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5, ncol = 3)

# Subset data (these parameters by default, modify thresholds to filter differently based on your data)
# QC filtering is performed based on the cell distribution depending on the above mentioned parameters.
# Check distribution of cells in the Violin plot and keep the main block of cells with a reasonable
# amount of features detected (usually more than 200), and low mitochondrial contamination. 
scdata <- subset(scdata, subset = nFeature_RNA>200 & nFeature_RNA<5000 & percent.mt<20)

# Normalization
############################################################
# We employ a global-scaling normalization method "LogNormalize" that normalizes the feature 
# expression measurements for each cell by the total expression, multiplies this by a scale factor 
# (10,000 by default), and log-transforms the result. Normalized values are stored in scdata[["RNA"]]@data.
scdata <- NormalizeData(scdata, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
############################################################
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset
# (i.e,  are highly expressed in some cells, and lowly expressed in others). Focusing on these genes 
# in downstream analysis helps to highlight biological signal in single-cell datasets.
# By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)

# Add sample information to each cell
scdata$sampletype <- colnames(scdata)
for (i in 1:length(sampletype)){
  scdata@meta.data[grep(paste0("-", i), colnames(scdata)), "sampletype"] <- sampletype[i]
}
#Some times you can pool sample types into conditions
#For example, if you have several treated and several control samples
#you can build 2 conditions, TREATED and CTL:
scdata$condition <- gsub("^.*?_", "", scdata$sampletype)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(scdata)
scdata <- ScaleData(scdata, verbose = TRUE, features = all.genes)
scdata <- RunPCA(scdata, verbose = TRUE)
DimPlot(scdata, reduction = "pca", dims = c(1,2))

# t-SNE and uMAP
scdata <- RunTSNE(scdata, reduction = "pca", dims = 1:20)
DimPlot(scdata, reduction = "tsne", dims = c(1,2))
scdata <- RunUMAP(scdata, reduction = "pca", dims = 1:20)
DimPlot(scdata, reduction = "umap", dims = c(1,2))

# Clustering
#In order to find clusters, we use a external function named kmeans
#This is a particular algorithm that finds clusters for a predefined number of clusters
#Define the number of clusters you expect to find:
numberofclusters <- 3
#You should perform clustering analysis with different number of clusters and decide which is the
#correct (real) number of clusters as the last number before the algorithm starts finding clusters
#that give no information (small group of some cells that are found as a cluster)

# The clustering analysis is perfomed differently every time you run the algorithm. By default we set
# the seed at 1, but try the analysis different times to ensure the clustering you performed is sensible.
set.seed(1)
scdata$kmeans <- kmeans(scdata@reductions[["pca"]]@cell.embeddings, centers = numberofclusters)$cluster
TSNEPlot(scdata, group.by = "kmeans") + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank()) + coord_flip() + scale_x_reverse()

ggsave(
  paste0("figs/tSNE_kmeans_clustering_", numberofclusters, ".pdf"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

# You can determine the factor by which you color/group the cells using group.by
#You can color/group cells by condition:
TSNEPlot(scdata, group.by = "condition", order = "CTL", 
         cols = c("#00BFC4", "#F8766D"), pt.size = 0.025) + theme(
           axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank()) + coord_flip() + scale_x_reverse()

ggsave(
  paste0("figs/tSNE_by_condition.pdf"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 5,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

# Or you can even split the plot by condition, to see if some cells only appear in one condition:
# Also you can use shape.by to plot different point shapes for each condition.
TSNEPlot(scdata, split.by = "condition", group.by = "kmeans",
         pt.size = 0.025) + theme(
           axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.title = element_blank(),
           axis.ticks = element_blank()) + coord_flip() + scale_x_reverse()

ggsave(
  paste0("figs/tSNE_kmeans_clustering_", numberofclusters, "_by_condition", ".pdf"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 5.5,
  height = 3,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)

# Visualization
#Feature plots are an interesting visualization strategy that allows representation of the expression
#of a single gene along all cells in the experiment. You can visualize genes of interest that can
#help you identify the phenotype of the different clusters you see.

#Read your txt containing the list of genes of interest:
features <- read.table("data/Table1.txt", header = F, stringsAsFactors = F)$V1


p <- FeaturePlot(scdata, reduction = "tsne",
                 features = features, 
                 #if there are too many features to plot at once,
                 #plot first the first half of features and then the rest
                 max.cutoff = 3,
                 cols = c("grey", "red"), 
                 combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend() + coord_flip() + scale_x_reverse()
}

cowplot::plot_grid(plotlist = p, ncol = 6)
#Determine the number of columns you want to distribute the plots in, depending on the total number
#of plots you want to make.

ggsave(
  "figs/feature_plots.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 7,
  height = 1.67,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)


#Another type of visualization are heatmaps. You can see gene expression on each cluster for a set
#genes of interest in the format of a heatmap
Idents(scdata) <- "kmeans" #To determine the kmeans clustering as the feature to group cells in the heatmap
Idents(scdata) <- factor(Idents(scdata), levels = 1:length(levels(Idents(scdata))))
DoHeatmap(subset(scdata, downsample = 1000), features = features[features %in% rownames(scdata)],
          size = 3, angle = 0) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")

#If you know the phenotype of the clusters, you can change the names of the Idents in order to
#show the cluster names in the heatmap
cluster.names <- c("Th0", "TFH", "TR1")
names(cluster.names) <- levels(scdata)
scdata <- RenameIdents(scdata, cluster.names)

DoHeatmap(subset(scdata, downsample = 1000), features = features[features %in% rownames(scdata)],
          size = 3, angle = 0) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")

ggsave(
  "figs/Heatmap.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 8,
  units = "in",
  dpi = 300,
  limitsize = TRUE
)


# Finding differentially expressed features (cluster biomarkers)
################################################################
# By default, it identifes positive and negative markers of a single cluster (specified in ident.1), 
# compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also
# test groups of clusters vs. each other, or against all cells.

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# By default, FindAllMarkers and FindMarkers function show only genes with a logfc threshold of 0.25 and
# min.pct of 0.1. If you want to recover more (but less significant) genes, you can modify these thresholds
Idents(scdata) <- "kmeans"
markers <- FindAllMarkers(scdata,
                          test.use = "negbinom", min.pct = 0.01)
write.table(markers, "output/Cluster_markers.txt", quote = F, sep = "\t")

#You can compare cluster or conditions 1 vs 1:
Idents(scdata) <- "condition"
COND <- FindMarkers(scdata, ident.1 = levels(Idents(scdata))[2], ident.2 = levels(Idents(scdata))[1],
                    verbose = T, test.use = "negbinom")

write.table(COND, "output/Condition_markers.txt", quote = F, sep = "\t")
