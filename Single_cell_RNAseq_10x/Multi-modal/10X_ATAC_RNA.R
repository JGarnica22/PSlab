# WNN analysis of 10x Multiome RNA+ATAC

# This script analyze multi-modal data consisting of gene expression and ATAC
# using  ‘weighted-nearest neighbor’ (WNN) analysis, an unsupervised framework
# to learn the relative utility of each data type in each cell, enabling an 
# integrative analysis.

# The workflow consists of three steps:
    # Independent preprocessing and dimensional reduction of each modality individually
    # Learning cell-specific modality ‘weights’, and constructing a WNN graph that
      # integrates the modalities
    # Downstream analysis (i.e. visualization, clustering, etc.) of the WNN graph

# https://satijalab.org/seurat/

# Load required packages:
cran.packages <- c("Seurat", "Signac", "tidyverse") #Best last version
for (i in cran.packages) {
  if(!require(i, character.only = T)) {
    install.packages(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
  library(i, character.only = T)
}

# Set here the species you are working with:
species <- "mouse" # say here either "mouse" or "human"
if (species == "mouse") {
  ensdb <- "EnsDb.Mmusculus.v79"
  genome <- "mm10" # indicate here the genome used for mapping!!
} else {
  ensdb <- "EnsDb.Hsapiens.v86"
  genome <- "hg38" # indicate here the genome used for mapping!!
}
bioc.packages <- c(ensdb)
for (i in bioc.packages) {
  if (!require(i, character.only = T)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
  library(i, character.only = T)
}


# First of all, we will load the output data from cellranger-arc
# We recommend downloading everything if possible, however,
# make sure to have at least these files:
    # Filtered feature barcode matrix (HDF5)
    # ATAC Per fragment information file (TSV.GZ)
    # ATAC Per fragment information index (TSV.GZ index)

# Then we will create a Seurat object based on the gene expression data, and 
# then add in the ATAC-seq data as a second assay.
setwd("C:/Users/Garnica/OneDrive - Universitat de Barcelona/PS")

# load 10x hdf5 file containing both RNA and ATAC data. 
data <- Read10X_h5("data/PS_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks

# Create Seurat object
scdata <- CreateSeuratObject(counts = rna_counts)
if (species == "mouse") {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^mt-")
} else {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
}

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb) #extract gene annotations from EnsDb
seqlevelsStyle(annotations) <- 'UCSC' #change to UCSC style (chr?)
genome(annotations) <- genome


# Load ATAC Per fragment information file (TSV.GZ)
frag.file <- "data/atac_fragments.tsv.gz"
scdata[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = genome,
  fragments = frag.file,
  min.cells = 10, # Subset accordingly!
  annotation = annotations
)


# perform basic QC based on the number of detected molecules for each modality
# as well as mitochondrial percentage.
VlnPlot(scdata, features = c("nCount_ATAC", "nCount_RNA","percent.mt"),
        ncol = 3, log = TRUE, pt.size = 0.5) + NoLegend()

# Subset based on QC plots!
scdata <- subset(
  x = scdata,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)


# We next perform pre-processing and dimensional reduction on both assays 
# independently, using standard approaches for RNA and ATAC-seq data.
# Remember to specify the assay before!

# RNA analysis
DefaultAssay(scdata) <- "RNA"
# Attention SCTtransform function normalize data, so input data must be non-normalized
scdata <- SCTransform(scdata, verbose = FALSE) %>%
          RunPCA() %>%
          RunTSNE(dims = 1:50, reduction.name = 'tsne.rna',
                  reduction.key = 'rnatsne_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with
# sequencing depth. This function also normalize data.
DefaultAssay(scdata) <- "ATAC"
scdata <- RunTFIDF(scdata)
scdata <- FindTopFeatures(scdata, min.cutoff = 'q0')
scdata <- RunSVD(scdata)
scdata <- RunTSNE(scdata, reduction = 'lsi', dims = 2:50,
                  reduction.name = "tsne.atac", reduction.key = "atactsne_")

# Next we calculate a WNN graph, representing a weighted combination of RNA and
# ATAC-seq modalities. We use this graph for TSNE visualization and clustering

scdata <- FindMultiModalNeighbors(scdata, reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:50, 2:50))
scdata <- RunTSNE(scdata, nn.name = "weighted.nn", reduction.name = "wnn.tsne",
                  reduction.key = "wnntsne_")
scdata <- FindClusters(scdata, graph.name = "wsnn", algorithm = 3,
                       verbose = FALSE)

# Visualize and explore graphs with clusters
plrna <- DimPlot(scdata, reduction = "tsne.rna") + ggtitle("RNA")
platac <- DimPlot(scdata, reduction = "tsne.atac") + ggtittle("ATAC")
plwnn <- DimPlot(scdata, reduction = "wnn.tsne") + ggtitle("RNA")
plrna + platac + plwnn & NoLegend() & theme(plot.title = element_text(hjust = 0.5))



# Plot aggregated signal
# https://satijalab.org/signac/articles/visualization.html

# Here is scripted a loop to represent multiple graphs for a list of genomic
# regions.
#Read your txt containing the list of genes or regions of interest:
features <- read.table("data/Table1.txt", header = F,
                       stringsAsFactors = F)$V1


# set how you want to group your data in the following plots:
Idents(scdata) <- "seurat_clusters" # could be sample type, condition... too
#This loop may take a while...
#Generate a pdf file for each feature.
for (f in features){
  pdf(paste0("figs/atac_combi_plots_", f, ".pdf"), width = 8, height = 12)
    # Compute the averaged frequency of sequenced DNA fragments for different
    # groups of cells within a given genomic region.  
    cov_plot <- CoveragePlot(
        object = scdata,
        region = f, #here you can use gene name or genomic region
        annotation = FALSE,
        peaks = FALSE)
    #Gene annotations within a given genomic region
    gene_plot <- AnnotationPlot(
      object = scdata,
      region = f)
    #Peak coordinates within a genomic region
    peak_plot <- PeakPlot(
      object = scdata,
      region = f)
    # Relationships between genomic positions: display an arc connecting two
    # linked positions, with the transparency of the arc line proportional to
    # a score associated with the link
    link_plot <- LinkPlot(
      object = scdata,
      region = f)
    # Inspect the frequency of sequenced fragments within a genomic region for 
    # individual cells, without aggregation in contrast to CoveragePlot()
    tile_plot <- TilePlot(
      object = scdata,
      region = f)
    # Visualize gene expression (RNA) data for each identity
    expr_plot <- ExpressionPlot(
      object = scdata,
      features = f,
      assay = "RNA"
    )
    # Combine previous plots into a single one
    CombineTracks(
      plotlist = list(cov_plot, tile_plot, peak_plot, gene_plot, link_plot),
      expression.plot = expr_plot,
      heights = c(10, 6, 1, 2, 3),
      widths = c(10, 1) )
  dev.off()
}



