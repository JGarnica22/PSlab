#This script is to analyze Smartseq2 scRNAseq data
#It was created for R 3.6.3 version (2021-01-12)
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
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(magrittr)
library(readxl)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(stringr)

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Terminal/SmartSeq2/Seurat")

# Specify parameters to be used along the script:
#Indicate species working with (typically mouse or human)
species <- "mouse"
#Indicate samples you sequenced separatedly (for example, you sorted your treated
#and control populations separatedly):
#It's better to indicate sample type as TREATMENT_CONDITION
#Load your sample-condition dataframe in order to link every cell to each condition:
samples <- read_xls(paste0(getwd(),"/data/SANTAMARIA_03.xls"), col_names = T,skip = 1)
names(samples)[10] <- "name"
samples <- samples[order(samples$name),] %>% as.data.frame()
rownames(samples) <- samples$name


# Read the data:
# To convert your data into a Seurat objects you should import a dataframe which each row names is a gene and each column
# represents a cell.

all <- read.table("data/Reads_all_samples.txt", sep = "\t", quote = "",
                  dec = ".", header = T, na.strings = T)
# Discard all rows with only 0 in all the cells
all2 <- all[rowSums(all)>0, ]
# Arrange cell names
names(all2)[671:length(names(all2))] <- sapply(strsplit(names(all2)[671:length(names(all2))], "_1", fixed = T), "[", 1)
names(all2) <- str_replace_all(names(all2), "[.]" , "-")



# if using STAR or other programs that use ensembl_id, convert them to gene_with biomaRt
en <- sapply(strsplit(rownames(all2), ".", fixed = T), "[", 1) %>% as.data.frame()
names(en) <- "ensembl_gene_id"
if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
}
BM <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
            mart = ensembl, verbose = T)
enbm <- merge(en, BM, by="ensembl_gene_id")
row.names(all2) <- make.names(enbm$external_gene_name, unique = T)

# We next use the count matrix to create a Seurat object. The object serves as a 
# container that contains both data (like the count matrix) and analysis (like PCA, 
# or clustering results) for a single-cell dataset.
    # min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. 
    # To reintroduce excluded features, create a new object with a lower cutoff.
    # min.features: Include cells where at least this many features are detected.

scdata <- CreateSeuratObject(counts = all2, min.cells = 3, min.features = 200)

# Add sample, condition and batch or project information to each cell
# Just need to add a new column to metadata
scdata@meta.data$sampletype <- NA
for (i in rownames(scdata@meta.data)){
  scdata@meta.data[i , "sampletype"] <- sub(".*?_", "", samples[ i ,"SAMPLE NAME"])
}
scdata@meta.data[is.na(scdata@meta.data$sampletype), "sampletype"] <- "Tet+_BDC"


# conditions
scdata@meta.data$condition <- NA
for (t in c("Tet[+]","Tet[-]")){
  scdata@meta.data[grep(t, scdata@meta.data$sampletype),"condition"] <- t
}


## Add TCR clonotype information
tracer <- read.csv(paste0(getwd(),"/data/cell_data.csv"))

tracer <- tracer %>% separate(A_productive, c("A_produc_TRV", "A_produc_CDR3", "A_produc_TRJ"), "_", remove = F)
tracer <- tracer %>% separate(B_productive, c("B_produc_TRV", "B_produc_CDR3", "B_produc_TRJ"), "_", remove = F)
scdata@meta.data$cell_name <- rownames(scdata@meta.data)
scdata@meta.data <- merge(scdata@meta.data, tracer, by = "cell_name", all = T )
rownames(scdata@meta.data) <- scdata@meta.data$cell_name
scdata@meta.data <- scdata@meta.data %>% dplyr::select(-cell_name) %>% filter(is.na(sampletype) == F)


# Perform quality control of Selection and filtration of cells based on QC metrics
# before integration as it needs to be done on raw RNA counts.

# QC and selecting cells for further analysis
#Criteria:
# * The number of unique genes detected in each cell (Low-quality cells or empty droplets will 
# often have very few genes and cell doublets or multiplets may exhibit an aberrantly high gene count
# * The percentage of reads that map to the mitochondrial genome (Low-quality/dying cells often
# exhibit extensive mitochondrial contamination

# Calculate mitochondrial QC metrics with the PercentageFeatureSet function:
if (species == "mouse") {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^mt.")
} else {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT.")
}

# Set sampletype as orig.ident to group cells in order to visualize the QC
Idents(scdata) <- scdata@meta.data$sampletype

# QC metrics can be found in the Seurat object metadata:
head(scdata@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5, ncol = 3)

# Subset data (these parameters by default, modify thresholds if you like to filter differently)
scdata <- subset(scdata, subset = nFeature_RNA>300 & nFeature_RNA<2000 & percent.mt<10)

# Run workflow for visualization and clustering normally
all.genes <- rownames(scdata)
scdata@meta.data$clonal_group <- as.factor(scdata@meta.data$clonal_group)
scdata <- NormalizeData(scdata, normalization.method = "LogNormalize", scale.factor = 10000)
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
scdata <- ScaleData(scdata, verbose = TRUE, features = all.genes)
scdata <- RunPCA(scdata, verbose = TRUE)
DimPlot(scdata, reduction = "pca", dims = c(1,2))

# t-SNE and uMAP
scdata <- RunTSNE(scdata, reduction = "pca", dims = 1:20)
DimPlot(scdata, reduction = "tsne", dims = c(1,2)) + DimPlot(scdata, group.by = "project", reduction = "tsne", dims = c(1,2))
DimPlot(scdata, group.by = "A_produc_TRV", reduction = "tsne", dims = c(1,2))
scdata <- RunUMAP(scdata, reduction = "pca", dims = 1:20)
DimPlot(scdata, reduction = "umap", dims = c(1,2)) + DimPlot(scdata, group.by = "project", reduction = "umap", dims = c(1,2))


## seurat clustering
scdata <- FindNeighbors(scdata, dims = 1:10)
scdata <- FindClusters(scdata, resolution = 0.5)
DimPlot(scdata, reduction = "tsne", dims = c(1,2))
TSNEPlot(scdata, group.by = "seurat_clusters", pt.size = 3) + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank()) + coord_flip() + scale_x_reverse()


# Clustering
#In order to find clusters, we use a external function named kmeans
#This is a particular algorithm that finds cluster for an already defined number of clusters
#Define the number of clusters you expect to find:
numberofclusters <- 3
#You should perform clustering analysis with different number of clusters and decide which is the
#correct (real) number of clusters as the last number before the algorithm starts finding clusters
#that give no information (small group of some cells that are found as a cluster)

set.seed(1)
scdata$kmeans <- kmeans(scdata@reductions[["pca"]]@cell.embeddings, centers = numberofclusters)$cluster
scdata@meta.data$kmeans <- as.factor(scdata@meta.data$kmeans)
TSNEPlot(scdata, group.by = "kmeans", pt.size = 3) + theme(
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
TSNEPlot(scdata, group.by = "sampletype", 
          pt.size = 3) + theme(
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
TSNEPlot(scdata, split.by = "kmeans", group.by = "sampletype", shape.by = "condition",
         pt.size = 3) + theme(
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



# Finding differentially expressed features (cluster biomarkers)
################################################################
# By default, it identifes positive and negative markers of a single cluster (specified in ident.1), 
# compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also
# test groups of clusters vs. each other, or against all cells.

# Find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(scdata) <- "kmeans"
DefaultAssay(scdata) <- "RNA"
markers <- FindAllMarkers(scdata,
                          test.use = "negbinom", min.pct = 0.05)
write.table(markers, "output/Cluster_markers.txt", quote = F, sep = "\t")

#You can compare cluster or conditions 1 vs 1:
Idents(scdata) <- "condition"
COND <- FindMarkers(scdata, ident.1 = levels(Idents(scdata))[2], ident.2 = levels(Idents(scdata))[1],
                    verbose = T, test.use = "negbinom",
                    logfc.threshold = 0, min.pct = 0.01)

write.table(COND, "output/Condition_markers.txt", quote = F, sep = "\t")



#### Compare with other files
prev <- readRDS(paste0(getwd(), "/previous/RDS/SANTAMARIA_01+03_clustered.rds"))

