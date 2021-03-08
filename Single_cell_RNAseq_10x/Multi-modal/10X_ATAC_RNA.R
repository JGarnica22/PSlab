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
cran.packages <- c("Seurat", "Signac", "tidyverse",
                   "patchwork", "devtools") #Best last version
for (i in cran.packages) {
  if(!require(i, character.only = T)) {
    install.packages(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
  library(i, character.only = T)
}

devtools.packages <- c("immunogenomics/presto") #Best last version
for (i in devtools.packages) {
  if(!require(i, character.only = T)) {
    devtools::install_github(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
  library(i, character.only = T)
}
devtools::install_github("immunogenomics/presto")
#############################################################################
## Set here the species you are working with:
#############################################################################
species <- "mouse" # say here either "mouse" or "human"
if (species == "mouse") {
  ensdb <- "EnsDb.Mmusculus.v79"
  bsgen <- "BSgenome.Mmusculus.UCSC.mm10"
  genome <- "mm10" # indicate here the genome used for mapping!!
  mt <- "^mt-"
} else {
  ensdb <- "EnsDb.Hsapiens.v86"
  bsgen <- "BSgenome.Hsapiens.UCSC.hg38"
  genome <- "hg38" # indicate here the genome used for mapping!!
  mt <- "^MT-"
}
bioc.packages <- c(ensdb, bsgen, "JASPAR2020", "TFBSTools", "motifmatchr",
                   "chromVAR")
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

###########################################################################
## Indicate here the folder where to find the cellranger-arc data
###########################################################################
setwd("C:/Users/Garnica/OneDrive - Universitat de Barcelona/PS")

# load 10x hdf5 file containing both RNA and ATAC data. 
data <- Read10X_h5("data/PS_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks

# Create Seurat object
scdata <- CreateSeuratObject(counts = rna_counts)
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = mt)


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
scdata <- SCTransform(scdata, verbose = FALSE) %>% #Perform NormalizeData, FindVariableFeatures,
                                                   # and ScaleData workflow
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
plwnn <- DimPlot(scdata, reduction = "wnn.tsne") + ggtitle("WNN")
plrna + platac + plwnn & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# You can peform sub-clustering on specifics clusters to find additional structures:
scdata <- FindSubCluster(scdata, cluster = 6, graph.name = "wsnn", algorithm = 3)

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



# Next we will examine the accessible regions of each cell to determine
# enriched motifs. We will use the chromVAR package from the Greenleaf lab.
# This calculates a per-cell accessibility score for known motifs, and adds 
# these scores as a third assay (chromvar) in the Seurat object.

# Scan the DNA sequence of each peak for the presence of each motif, and create
# a Motif object
DefaultAssay(scdata) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(scdata), pwm = pwm_set,
                                  genome = genome, use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
scdata <- SetAssayData(scdata, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes!!
scdata <- RunChromVAR(
  object = scdata,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Then, we explore the multimodal dataset to identify key regulators of each cell state.
# Paired data allow identification of transcription factors (TFs) 
# that satisfy multiple criteria, helping to narrow down the list of putative 
# regulators to the most likely candidates. We aim to identify TFs whose expression
# is enriched in multiple cell types in the RNA measurements, but also have enriched
# accessibility for their motifs in the ATAC measurements.

# We can then represent features plots of gene expression and accessibility enrichment
# of some TFs.
# First, convert from motif name to motif ID or vice versa to get motifs for a certain TF.

feat_TF <- features[] #only list TF!

pdf("figs/expression_motifs_plots_.pdf", width = 12, height = 8)
  for (f in feat_TF){
    motif.name <- ConvertMotifID(scdata, name = f) 
    gene_plot <- FeaturePlot(scdata, features = paste0("sct_", f),
                             reduction = 'wnn.tsne')
    motif_plot <- FeaturePlot(scdata, features = motif.name,
                              min.cutoff = 0, cols = c("lightgrey", "darkred"),
                              reduction = 'wnn.tsne')
    gene_plot | motif_plot
  }
dev.off()


# If we want to quantify this relationship and search across all cell types, 
# we will use the presto package to perform fast differential expression. 
# We run two tests: one using gene expression data, and the other using chromVAR 
# motif accessibilities. presto calculates a p-value based on the Wilcox rank sum test,
# which is also the default test in Seurat, and we restrict our search to TFs that
# return significant results in both tests.

# presto also calculates an “AUC” statistic, which reflects the power of each gene
# (or motif) to serve as a marker of cell type. A maximum AUC value of 1 indicates
# a perfect marker. Since the AUC statistic is on the same scale for both genes
# and motifs, we take the average of the AUC values from the two tests, and use
# this to rank TFs for each cell type:
  
markers_rna <- wilcoxauc(scdata,
                         group_by = 'seurat_clusters', #group by samples or condition!
                         assay = 'data',
                         seurat_assay = 'SCT')
markers_motifs <- wilcoxauc(scdata,
                            group_by = 'seurat_clusters', #group by samples or condition!
                            assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(pbmc, id = motif.names)

