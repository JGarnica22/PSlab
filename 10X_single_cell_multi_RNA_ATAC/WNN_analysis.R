# WNN analysis of 10x Multiome RNA+ATAC

# This script analyze multi-modal data consisting of gene expression and ATAC
# using  'weighted-nearest neighbor' (WNN) analysis, an unsupervised framework
# to learn the relative utility of each data type in each cell, enabling an 
# scdataegrative analysis.

# The workflow consists of three steps:
# Independent preprocessing and dimensional reduction of each modality individually
# Learning cell-specific modality 'weights', and constructing a WNN graph that
# scdataegrates the modalities
# Downstream analysis (i.e. visualization, clustering, etc.) of the WNN graph

# https://satijalab.org/seurat/

.libPaths("C:/Users/Garnica/OneDrive - Universitat de Barcelona/Documentos/R/win-library/4.0/")
# Load required packages:
cran.packages <- c("Seurat", "Signac", "tidyverse",
                   "patchwork", "devtools", "BioVenn") #Best last version
for (i in cran.packages) {
  # if(!require(i, character.only = T)) {
  #   install.packages(i)
  #   prscdata(paste(i,"just installed"))
  # } else {
  #   prscdata(paste(i,"was already installed"))
  # }
  library(i, character.only = T)
}

# devtools.packages <- c("immunogenomics/presto") #Best last version
# for (i in devtools.packages) {
#   # if(!require(i, character.only = T)) {
#   #   devtools::install_github(i)
#   #   prscdata(paste(i,"just installed"))
#   # } else {
#   #   prscdata(paste(i,"was already installed"))
#   # }
#   library(i, character.only = T)
# }
# devtools::install_github("immunogenomics/presto")
library(presto)
#############################################################################
## Set here the species you are working with:
#############################################################################
species <- "mouse" # say here either "mouse" or "human"
if (species == "mouse") {
  ensdb <- "EnsDb.Mmusculus.v79"
  bsgen <- "BSgenome.Mmusculus.UCSC.mm10"
  genome <- "mm10" # indicate here the genome used for mapping!!
  mt <- "^mt-"
  blacklist <- "blacklist_mm10"
  txid <- 10090
  
} else {
  ensdb <- "EnsDb.Hsapiens.v86"
  bsgen <- "BSgenome.Hsapiens.UCSC.hg38"
  genome <- "hg38" # indicate here the genome used for mapping!!
  mt <- "^MT-"
  blacklist <- "blacklist_hg38_unified"
  txid <- 9609
  
}
bioc.packages <- c(ensdb, bsgen, "JASPAR2020", "TFBSTools", "motifmatchr",
                   "chromVAR", 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                   'limma', 'S4Vectors', 'SingleCellExperiment',
                   'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                   "GenomicRanges", "Gviz", "trackViewer", "rtracklayer", 
                   "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db")
for (i in bioc.packages) {
  # if (!require(i, character.only = T)) {
  #   BiocManager::install(i)
  #   print(paste(i,"just installed"))
  # } else {
  #   print(paste(i,"was already installed"))
  # }
  library(i, character.only = T)
}
if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
}
titol <- function(x) {plot_annotation(title = x,
                                      theme = theme(plot.title = element_text(size = 16,
                                                                              hjust = 0.5, face = "bold")))}


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
setwd("~/Bioinformatics/inkt_multiome_santamaria_21_23/analysis_21")


# load 10x hdf5 file containing both RNA and ATAC data. 
data <- Read10X_h5("data/aggr_multiome/filtered_feature_bc_matrix.h5")


## !!!! Merge objects from different samples???
# extract RNA and ATAC data
rna_counts <- data$`Gene Expression`
atac_counts <- data$Peaks

# Create Seurat object
scdata <- CreateSeuratObject(counts = rna_counts)
# Create QC metrics for mithocondrial content
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = mt)


# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
library(GenomeInfoDb)
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = get(ensdb)) #extract gene annotations from EnsDb
seqlevelsStyle(annotations) <- 'UCSC' #change to UCSC style (chr?)
genome(annotations) <- genome


# Load ATAC Per fragment information file (TSV.GZ)
frag.file <- "data/aggr_multiome/atac_fragments.tsv.gz"
scdata[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = genome,
  fragments = frag.file,
  min.cells = 3, # Subset accordingly!
  min.features = 10,
  annotation = annotations
)


# scdata[["ATAC"]] <- CreateSeuratObject(
#   counts = scdata[["ATAC"]],
#   assay = "peaks",
#   meta.data = scdata@meta.data
# )
# library(trackViewer)
# pea <- importScore("data/treated/atac_peaks.bed", 
#                    format = "BED",
#                    ignore.strand = T)
# scdata[["peaks"]] <- CreateChromatinAssay(
#   counts = gr_obj,
#   fragments = frag.file,
#   annotation = annotation
# )
# 
#gr_obj =  import("data/aggr_multiome/atac_peaks.bed")


### Alternatively you can also call peaks using mac2:#####################
# The set of peaks identified using Cellranger often merges distinct peaks 
# that are close together. This can create a problem for certain analyses, 
# particularly motif enrichment analysis and peak-to-gene linkage. To identify 
# a more accurate set of peaks, we can call peaks using MACS2 with the CallPeaks()
# function. Here we call peaks on all cells together, but we could identify peaks 
# for each group of cells separately by setting the group.by parameter, and this 
# can help identify peaks specific to rare cell populations

# # call peaks using MACS2
# peaks <- CallPeaks(scdata, 
#                    macs2.path = "/home/josep/software/anaconda3/envs/macs2")
# # remove peaks on nonstandard chromosomes and in genomic blacklist regions
# peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
# peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)
# 
# # quantify counts in each peak
# macs2_counts <- FeatureMatrix(
#   fragments = Fragments(scdata),
#   features = peaks,
#   cells = colnames(scdata))
# 
# create a new assay using the MACS2 peak set and add it to the Seurat object
# scdata[["peaks"]] <- CreateChromatinAssay(
#   counts =gr_obj,
#   fragments = fragpath,
#   annotation = annotation
# )


scdata$sampletype <- colnames(scdata)
scdata@meta.data[grep("-1", colnames(scdata)), "sampletype"] <- "CTL"
scdata@meta.data[grep("-2", colnames(scdata)), "sampletype"] <- "NP"

# Create QC metrics for ATAC data
remove(data)
remove(rna_counts)
remove(atac_counts)

library(GenomicRanges)
library(GenomicAlignments)
DefaultAssay(scdata) <- "ATAC"
scdata <- NucleosomeSignal(scdata)
scdata <- TSSEnrichment(scdata)

# scdata$pct_reads_in_peaks <- scdata$peak_region_fragments / 
#                              scdata$passed_filters * 100
# scdata$blacklist_ratio <- scdata$blacklist_region_fragments /
#                           scdata$peak_region_fragments
# scdata$blacklist_fraction <- FractionCountsInRegion(
#                                                     object = scdata,
#                                                     assay = 'ATAC',
#                                                     regions = blacklist)
# DefaultAssay(scdata) <- "peaks"
# scdata <- NucleosomeSignal(scdata)
# scdata <- TSSEnrichment(scdata)

# perform basic QC based on the number of detected molecules for each modality
# as well as mitochondrial percentage.

## LABEL DATA



preqc <- VlnPlot(scdata, group.by="sampletype", features = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                                             "nCount_ATAC", "nFeature_ATAC",
                                                             "TSS.enrichment", "nucleosome_signal"),
                 # "nCount_peaks", "nFeature_peaks",
                 # "pct_reads_in_peaks", "peak_region_fragments",
                 # "blacklist_ratio", "blacklist_fraction"),
                 ncol = 4, log = T, pt.size = 0.5) + NoLegend()

# Subset based on QC plots!!
scdata <- subset(
  x = scdata,
  subset = nCount_ATAC < 9.9e4 &
    nCount_ATAC > 7e2 &
    nCount_RNA < 10000 &
    nCount_RNA > 500 &
    percent.mt < 12 &
    nucleosome_signal < 2.8 &
    TSS.enrichment > 3
  #blacklist_ratio < 0.05
)

# scdata$high.tss <- ifelse(scdata$TSS.enrichment > 2, 'High', 'Low')
# TSSPlot(scdata, group.by = 'high.tss') + NoLegend()
# scdata$nucleosome_group <- ifelse(scdata$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# FragmentHistogram(object = scdata, group.by = 'nucleosome_group')

# We next perform pre-processing and dimensional reduction on both assays 
# independently, using standard approaches for RNA and ATAC-seq data.
# Remember to specify the assay before!

# RNA analysis
DefaultAssay(scdata) <- "RNA"
# Attention SCTtransform function normalize data, so input data must be non-normalized
scdata <- SCTransform(scdata, verbose = FALSE) %>% #Perform NormalizeData, FindVariableFeatures,
  # and ScaleData workflow, stored in "SCT" assay
  RunPCA() %>%
  RunTSNE(dims = 1:50, reduction.name = 'tsne.rna',
          reduction.key = 'rnatsne_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with
# sequencing depth. This function also normalize data.
DefaultAssay(scdata) <- "ATAC" #or use "peaks here
scdata <- RunTFIDF(scdata)
# Signac performs term frequency-inverse document frequency (TF-IDF) normalization.
# This is a two-step normalization procedure, that both normalizes across cells to
# correct for differences in cellular sequencing depth, and across peaks to give 
# higher values to more rare peaks
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

## Load original data to compare

data_origi <- Read10X(data.dir = "data/iNKT_original")
scdata_origi <- CreateSeuratObject(counts = data_origi, min.cells = 3, min.features = 200)
scdata_origi[["percent.mt"]] <- PercentageFeatureSet(scdata_origi, pattern = "^mt-")


VlnPlot(scdata_origi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.5, ncol = 3) + titol("Original experiment PRE-SUBSET")

scdata_origi <- subset(scdata_origi, subset = nFeature_RNA>1000 &
                         nFeature_RNA<4000 & percent.mt<20)
scdata_origi[["RBC"]] <- PercentageFeatureSet(scdata_origi, pattern = "Hb")
scdata_origi[["CD19"]] <- PercentageFeatureSet(scdata_origi, pattern = "Cd19")
scdata_origi[["CD20"]] <- PercentageFeatureSet(scdata_origi, pattern = "Ms4a1")
scdata_origi <- subset(scdata_origi, subset = RBC < 5 & CD19 < 0.005 & CD20 < 0.005)

scdata_origi$sampletype <- colnames(scdata_origi)
scdata_origi@meta.data[grep("-1", colnames(scdata_origi)), "sampletype"] <- "CTL_original"
scdata_origi@meta.data[grep("-2", colnames(scdata_origi)), "sampletype"] <- "CTL_original"
scdata_origi@meta.data[grep("-3", colnames(scdata_origi)), "sampletype"] <- "Cd1d_original"
scdata_origi@meta.data[grep("-4", colnames(scdata_origi)), "sampletype"] <- "Cd1d_original"

# scdata_origi <- SCTransform(scdata_origi, verbose = FALSE) %>% 
#   #Perform NormalizeData, FindVariableFeatures,
#   # and ScaleData workflow, stored in "SCT" assay
#                 RunPCA() %>%
#                 RunTSNE(reduction = "pca", dims = 1:30)
scdata_origi <- NormalizeData(scdata_origi, normalization.method = "LogNormalize",
                              scale.factor = 10000)
scdata_origi <- FindVariableFeatures(scdata_origi, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scdata_origi)
# Run the standard workflow for visualization and clustering
scdata_origi <- ScaleData(scdata_origi, verbose = TRUE, features = all.genes)
scdata_origi <- RunPCA(scdata_origi, verbose = TRUE)
scdata_origi <- RunTSNE(scdata_origi, reduction = "pca", dims = 1:30)

set.seed(2)
scdata_origi$kmeans <- kmeans(scdata_origi@reductions[["pca"]]@cell.embeddings, centers = 5)$cluster
TSNEPlot(scdata_origi, group.by = "kmeans", split.by="sampletype") + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank()) + coord_flip() + scale_x_reverse()

scdata_origi$celltype <- scdata_origi$kmeans
scdata_origi@meta.data[grep("1", scdata_origi$kmeans), "celltype"] <- "Cluster3"
scdata_origi@meta.data[grep("2", scdata_origi$kmeans), "celltype"] <- "Cluster1"
scdata_origi@meta.data[grep("3", scdata_origi$kmeans), "celltype"] <- "Cluster5"
scdata_origi@meta.data[grep("4", scdata_origi$kmeans), "celltype"] <- "Cluster2"
scdata_origi@meta.data[grep("5", scdata_origi$kmeans), "celltype"] <- "Cluster4"

cols <- c("#A3A500", "#00B0F6", "#E76BF3", "#00BF7D", "#F8766D")

##  PREDICTION WITH ORIGINAL DATA
DefaultAssay(scdata) <- "RNA"
anchors <- FindTransferAnchors(reference = scdata_origi , query = scdata)
predictions <- TransferData(anchorset = anchors, refdata = scdata_origi$celltype, 
                            dims = 1:30)
scdata <- AddMetaData(scdata, metadata = predictions)

scdata$predicted.id <- factor(x = scdata$predicted.id, 
                              levels = sort(unique(scdata$predicted.id)))
remove(data_origi)
# Visualize and explore graphs with clusters
plrna_1 <- DimPlot(scdata[, scdata$sampletype == "CTL"], reduction = "tsne.rna", group.by = "predicted.id",
                   split.by="sampletype",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8))+ NoLegend()
plrna_2 <- DimPlot(scdata[, scdata$sampletype == "NP"], reduction = "tsne.rna", group.by = "predicted.id",
                   split.by="sampletype",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8)) + NoLegend()
plrna_3 <- DimPlot(scdata, reduction = "tsne.rna", group.by = "predicted.id",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8))

platac_1 <- DimPlot(scdata[, scdata$sampletype == "CTL"], reduction = "tsne.atac", group.by = "predicted.id",
                    split.by="sampletype",
                    cols = cols
) + NoAxes() + theme(title = element_text(size = 8)) + NoLegend()
platac_2 <- DimPlot(scdata[, scdata$sampletype == "NP"], reduction = "tsne.atac", group.by = "predicted.id",
                    split.by="sampletype",
                    cols = cols
) + NoAxes() + theme(title = element_text(size = 8)) + NoLegend()
platac_3 <- DimPlot(scdata, reduction = "tsne.atac", group.by = "predicted.id",
                    cols = cols
) + NoAxes() + theme(title = element_text(size = 8))
plwnn_1 <- DimPlot(scdata[, scdata$sampletype == "CTL"], reduction = "wnn.tsne", group.by = "predicted.id",
                   split.by="sampletype",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8)) + NoLegend()
plwnn_2 <- DimPlot(scdata[, scdata$sampletype == "NP"], reduction = "wnn.tsne", group.by = "predicted.id",
                   split.by="sampletype",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8)) + NoLegend()
plwnn_3 <- DimPlot(scdata, reduction = "wnn.tsne", group.by = "predicted.id",
                   cols = cols
) + NoAxes() + theme(title = element_text(size = 8))
library(ggpubr)
plrna <- ggarrange(plrna_1, plrna_2, plrna_3, ncol=3, widths = c(1,1,1.4)) %>% 
  annotate_figure(top = text_grob("RNA", color = "grey",
                                  face = "bold", size = 16))
platac <- ggarrange(platac_1, platac_2, platac_3, ncol=3, widths = c(1,1,1.4)) %>% 
  annotate_figure(top = text_grob("ATAC", color = "grey",
                                  face = "bold", size = 16))
plwnn <- ggarrange(plwnn_1, plwnn_2, plwnn_3, ncol=3, widths = c(1,1,1.4)) %>% 
  annotate_figure(top = text_grob("WNN", color = "grey",
                                  face = "bold", size = 16))




# You can peform sub-clustering on specifics clusters to find additional structures:
# scdata <- FindSubCluster(scdata, cluster = 6, graph.name = "wsnn", algorithm = 3)

## PDF only for tSNE
pdf("figs/tSNE_plots_2.pdf", width = 16, height = 5)
plrna
platac
plwnn
dev.off()


# Plot aggregated signal
# https://satijalab.org/signac/articles/visualization.html

# Here is scripted a loop to represent multiple graphs for a list of genomic
# regions.
#Read your txt containing the list of genes or regions of scdataerest:
features <- read.table("data/Table 2.txt", header = F,
                       stringsAsFactors = F)$V1


# set how you want to group your data in the following plots:
Idents(scdata) <- "predicted.id" # could be sample type, condition... too
#This loop may take a while...
#Generate a pdf file for each feature.
DefaultAssay(scdata) <- "ATAC"
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
    assay = "SCT"
  )
  # Combine previous plots scdatao a single one
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
library(TFBSTools)
DefaultAssay(scdata) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = txid,
                                                    all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(scdata), pwm = pwm_set,
                                  genome = get(bsgen), use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
scdata <- SetAssayData(scdata, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes!!
scdata <- RunChromVAR(
  object = scdata,
  genome = get(bsgen)
)


#Alternatively, find overepresetned motifs by:
# https://satijalab.org/signac/articles/motif_vignette.html

da_peaks <- FindAllMarkers(
  object = scdata,
  features = features,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC') #or peaks




# Then, we explore the multimodal dataset to identify key regulators of each cell state.
# Paired data allow identification of transcription factors (TFs) 
# that satisfy multiple criteria, helping to narrow down the list of putative 
# regulators to the most likely candidates. We aim to identify TFs whose expression
# is enriched in multiple cell types in the RNA measurements, but also have enriched
# accessibility for their motifs in the ATAC measurements.

# We can then represent features plots of gene expression and accessibility enrichment
# of some TFs.
# First, convert from motif name to motif ID or vice versa to get motifs for a certain TF.

feat_TF <- features #only list TF!

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

# presto also calculates an "AUC" statistic, which reflects the power of each gene
# (or motif) to serve as a marker of cell type. A maximum AUC value of 1 indicates
# a perfect marker. Since the AUC statistic is on the same scale for both genes
# and motifs, we take the average of the AUC values from the two tests, and use
# this to rank TFs for each cell type:


markers_rna <- wilcoxauc(scdata,
                         group_by = 'predicted.id', #group by samples or condition!
                         assay = 'data',
                         seurat_assay = 'SCT')
markers_motifs <- wilcoxauc(scdata,
                            group_by = 'predicted.id', #group by samples or condition!
                            assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(scdata, id = motif.names)



#############################################################################
## LINK PEAKS TO GENES ##
# For each gene, we can find the set of peaks that may regulate the gene by by 
# computing the correlation between gene expression and accessibility at nearby peaks,
# and correcting for bias due to GC content, overall accessibility, and peak size.

# Running this step on the whole genome can be time consuming, so you can use a
# subset of genes as an example. 

DefaultAssay(scdata) <- "ATAC"  #or "peaks"

# first compute the GC content for each peak
scdata <- RegionStats(scdata, genome = get(bsgen))

# link peaks to genes
scdata <- LinkPeaks(
  object = scdata,
  peak.assay = "ATAC", #or "peaks"
  expression.assay = "SCT",
  genes.use = features) #or remove this option if you want to run over all the genome.

# Visualize with coverageplot as done before.
pdf("figs/clustering_and_coverage_WNN.pdf", width = 12, height = 10)
print(plrna)
print(platac)
print(plwnn)
for (f in features[-40]){
  print(CoveragePlot(
    object = scdata,
    region = f,
    features = f,
    # assay= "ATAC",
    peaks = F,
    expression.assay = "SCT", tile = T,links = T,
    
    #idents = "seurat_clusters", #or sample type, conditions, etc
    extend.upstream = 10000, 
    extend.downstream = 10000))
}
dev.off()

CoverageBrowser(
  object = scdata,
  region = f,
  # assay= "ATAC",
)


## PERFORM THE CLASSICAL ANALYSIS FOR RNA AND ATAC
DefaultAssay(scdata) <- "RNA"
heatmap_combi <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                           features = features, slot = "data",
                           size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("Gene expression comparison")


DefaultAssay(scdata) <- "RNA"
p <- FeaturePlot(scdata, reduction = "wnn.tsne",
                 features = features, 
                 max.cutoff = 3,
                 cols = c("grey", "red"), 
                 combine = FALSE)

for(o in 1:length(p)) {
  p[[o]] <- p[[o]] + NoAxes() + NoLegend()
  
}


featur <- cowplot::plot_grid(plotlist = p, ncol = 7) +
  titol("iNKT markers list. WNN clustering")

pdf("figs/RNA_analysis.pdf", width = 12, height = 8)
heatmap_combi
featur
dev.off()

heat <- function(sample, features, titol) {
  DoHeatmap(subset(scdata[, scdata$sampletype  %in% sample],
                   downsample = 400),
            group.by = "predicted.id",
            features = features, slot = "data",
            size = 3, angle = 90) + NoLegend() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + titol(titol)
}

## version of heatmap
DefaultAssay(scdata) <- "RNA"
heat_rna <- ggarrange(heat("CTL", features, "CTL"),
                heat("NP", features, "NP") + 
                  theme(axis.text.y=element_blank()),
                heat(c("CTL","NP"), features, "ALL") + 
                  theme(axis.text.y=element_blank()),
                ncol = 3) %>% 
  annotate_figure(top = text_grob("RNA expression comparison",
                                  color = "grey",
                                  face = "bold", size = 16))
heat_atacrna <- ggarrange(heat("CTL", clo$gene_name, "CTL"),
                      heat("NP", clo$gene_name, "NP") + 
                        theme(axis.text.y=element_blank()),
                      heat(c("CTL","NP"), clo$gene_name, "ALL") + 
                        theme(axis.text.y=element_blank()),
                      ncol = 3, widths = c(1.2,1,1)) %>% 
  annotate_figure(top = text_grob("Expression of closest gene to top 50 OCR",
                                  color = "grey",
                                  face = "bold", size = 16))
h1 <- heat("CTL", features, "CTL RNA") 
h3 <- heat("NP", features, "NP RNA") + theme(axis.text.y=element_blank())
h5 <- heat(c("CTL","NP"), features, "ALL RNA") + 
  theme(axis.text.y=element_blank())
DefaultAssay(scdata) <- "ATAC"
heat_atac <- ggarrange(heat("CTL", top$gene, "CTL"),
                      heat("NP", top$gene, "NP") + 
                        theme(axis.text.y=element_blank()),
                      heat(c("CTL","NP"), top$gene, "ALL") + 
                        theme(axis.text.y=element_blank()),
                      ncol = 3, widths = c(1.4,1,1)) %>% 
  annotate_figure(top = text_grob("50 most differintally open OCR",
                                  color = "grey",
                                  face = "bold", size = 16))

closest_features <- data.frame(matrix(nrow = length(features), 
                                      ncol = length(clo_fea)))
names(closest_features) <- names(clo_fea)
for (i in 1:length(features)){
  if (features[i] != "Izumo1r"){
    closest_features[i,] <- clo_fea %>% filter(gene_name == features[i]) %>% 
      filter(distance == min(distance))
  }
}
closest_features <- closest_features[complete.cases(closest_features),]


h2 <- heat("CTL", closest_features$query_region, "CTL ORC") + 
  theme(axis.text.y=element_blank())
h4 <- heat("NP", closest_features$query_region, "NP ORC") + 
        theme(axis.text.y=element_blank())
h6 <- heat(c("CTL","NP"), closest_features$query_region, "ALL ORC")+ 
  theme(axis.text.y=element_blank())

heat_com <- ggarrange(h1,h2,h3,h4,h5,h6,
                      ncol = 6, widths = c(1.1,1,1,1,1)) %>% 
  annotate_figure(top = text_grob("Expressión and closest OCR accessability",
                                  color = "grey",
                                  face = "bold", size = 16))

tbl <- table(scdata$predicted.id, scdata$sampletype) %>% addmargins()

library(gridExtra)
pdf("figs/Heatmaps_2.pdf", width = 20, height = 10)
grid.table(tbl)
heat_rna
heat_atac
heat_atacrna
heat_com
dev.off()


DefaultAssay(scdata) <- "RNA"
Idents(scdata) <- "predicted.id"
markers_rna <- FindAllMarkers(scdata, test.use = "negbinom")
library(writexl)
writexl::write_xlsx(markers_rna, "out/markers_rna.xlsx")
markers_rna <- readxl::read_xlsx("out/markers_rna.xlsx") %>% 
  as.data.frame()
cl_1 <- markers_rna %>% dplyr::filter(cluster == "Cluster1")
cl_2 <- markers_rna %>% dplyr::filter(cluster == "Cluster2") 
cl_3 <- markers_rna %>% dplyr::filter(cluster == "Cluster3") 
cl_4 <- markers_rna %>% dplyr::filter(cluster == "Cluster4") 
cl_5 <- markers_rna %>% dplyr::filter(cluster == "Cluster5")
writeClipboard(cl_5$gene)





## ATAC
DefaultAssay(scdata) <- "ATAC"
markers_atac <- FindAllMarkers(scdata, test.use = "negbinom")
topc1 <-  subset(markers_atac, markers_atac$cluster == "Cluster5" & markers_atac$p_val_adj<=0.001)
top <- head(topc1[order(topc1$avg_log2FC, decreasing = T), ], 50)
clo <- ClosestFeature(scdata, regions = top$gene)
writexl::write_xlsx(markers_atac, "out/markers_atac.xlsx")
markers_atac <- readxl::read_xlsx("out/markers_atac.xlsx") %>% 
  as.data.frame()

cl_1 <- markers_atac %>% dplyr::filter(cluster == "Cluster1")
cl_2 <- markers_atac %>% dplyr::filter(cluster == "Cluster2") 
cl_3 <- markers_atac %>% dplyr::filter(cluster == "Cluster3") 
cl_4 <- markers_atac %>% dplyr::filter(cluster == "Cluster4") 
cl_5 <- markers_atac %>% dplyr::filter(cluster == "Cluster5")
writeClipboard(cl_5$gene)

p <- FeaturePlot(scdata, reduction = "wnn.tsne",
                 features = top$gene[1:25], 
                 max.cutoff = 3,
                 cols = c("grey", "red"), 
                 combine = FALSE)

for(o in 1:length(p)) {
  p[[o]] <- p[[o]] + NoAxes() + NoLegend()
}

featur_atac <- cowplot::plot_grid(plotlist = p, ncol = 7) +
  titol("iNKT markers list. WNN clustering")

p <- FeaturePlot(scdata, reduction = "wnn.tsne",
                 features = top$gene[26:50], 
                 max.cutoff = 3,
                 cols = c("grey", "red"), 
                 combine = FALSE)

for(o in 1:length(p)) {
  p[[o]] <- p[[o]] + NoAxes() + NoLegend()
}

featur_atac_2 <- cowplot::plot_grid(plotlist = p, ncol = 7) +
  titol("iNKT markers list. WNN clustering")

heatmap_combi_atac <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                                features = top$gene, slot = "data",
                                size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("50 most differenttialy opened sites in Cluster 5")

heatmap_combi_atac <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                                features = top$gene, slot = "data",
                                size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("50 most differentially opened sites in Cluster 5")
DefaultAssay(scdata) <- "RNA"
heatmap_combi_atac_clo <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                                    features = clo$gene_name, slot = "data",
                                    size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("Closest gene to top 50 sites in Cluster 5")

pdf("figs/ATAC_analysis.pdf", width = 12, height = 8)
heatmap_combi_atac
heatmap_combi_atac_clo
dev.off()


## Trajectories with monocle
library(SeuratWrappers)
library(monocle3)
library(Matrix)
DefaultAssay(scdata) <- "ATAC"

# iNKTR1 <- scdata[,  scdata$predicted.id == "Cluster5(iNKTR1)"]
inktr1.cds <- as.cell_data_set(scdata)
inktr1.cds <- preprocess_cds(inktr1.cds, method = "PCA")
inktr1.cds <- reduce_dimension(inktr1.cds, reduction_method = "UMAP")
inktr1.cds <- cluster_cells(cds = inktr1.cds, reduction_method = "UMAP")
inktr1.cds <- learn_graph(inktr1.cds, use_partition = TRUE)
DimPlot(scdata, reduction = "tsne.atac", group.by = "predicted.id") +
  ggtitle("ATAC-based clustering") + NoAxes()
inktr1.cds <- order_cells(inktr1.cds, reduction_method = "UMAP")


plot_cells(
  cds = inktr1.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

scdata <- Footprint(
  object = scdata,
  motif.name = features,
  genome = get(bsgen)
)

p2 <- PlotFootprint(bone, features = features)



DefaultAssay(scdata) <- "RNA"
heatmap_s1 <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                        features = features, slot = "data",
                        size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("Gene expression comparison")

DefaultAssay(scdata) <- "ATAC"
closest <- ClosestFeature(scdata, regions = rownames(scdata[["ATAC"]]))
clo_fea <- closest %>% dplyr::filter(gene_name %in% features) %>% 
  mutate(com = paste0(gene_name,"_",query_region))

heatmap_s1 <- DoHeatmap(subset(scdata, downsample = 1000), group.by = "predicted.id",
                        features = clo_fea$query_region, slot = "data", 
                        size = 3, angle = 55) + NoLegend() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  titol("Gene expression comparison")

library(ggpubr)
DefaultAssay(scdata) <- "ATAC"

for (x in features){
  pdf(paste0("figs/OCR_table1_", x, ".pdf"), width = 24, height = 12)
  cl <- clo_fea %>% dplyr::filter(gene_name == x)
  if (nrow(cl) > 1){
    pl <- ggarrange(heat("CTL", cl$query_region, "CTL"),
                    heat("NP", cl$query_region, "NP"),
                    heat(c("CTL","NP"), cl$query_region, "ALL"),
                    ncol = 3)
    pl <- annotate_figure(pl,
                          top = text_grob(paste0("OCR regions close to ", x),
                                          color = "grey",
                                          face = "bold", size = 16))
    print(pl)
    
  }  
  cl <- clo_fea %>% dplyr::filter(gene_name == x & distance <= 2000)
  if (nrow(cl) > 0){
    pl <- ggarrange(heat("CTL", cl$query_region, "CTL"),
                    heat("NP", cl$query_region, "NP"),
                    heat(c("CTL","NP"), cl$query_region, "ALL"),
                    ncol = 3)
    pl <- annotate_figure(pl,
                          top = text_grob(paste0("< 2000 kb far OCR regions ", x),
                                          color = "grey",
                                          face = "bold", size = 16))
    print(pl)
  }
  cl <- clo_fea %>% dplyr::filter(gene_name == x & distance > 2000)
  if (nrow(cl) > 0){
    pl <- ggarrange(heat("CTL", cl$query_region, "CTL"),
                    heat("NP", cl$query_region, "NP"),
                    heat(c("CTL","NP"), cl$query_region, "ALL"),
                    ncol = 3)
    pl <- annotate_figure(pl,
                          top = text_grob(paste0("> 2000 kb far OCR regions ", x),
                                          color = "grey",
                                          face = "bold", size = 16))
    print(pl)
    
  }
  dev.off()
}




Idents(scdata) <- "predicted.id"
DefaultAssay(scdata) <- "ATAC"
CoveragePlot(
  object = scdata,
  region = f,
  # features = f,
  # assay= "ATAC",
  peaks = F,
  big.wig = paste0("data/aggr_multiome/", bw[w]),
  bigwig.type = "coverage",
  extend.upstream = 10000, 
  extend.downstream = 10000)
DefaultAssay(scdata) <- "ATAC"
cells <- data.frame(rownames(scdata@meta.data),scdata@meta.data$predicted.id,
                    scdata@meta.data$sampletype)
cells2 <- cells %>% separate(rownames.scdata.meta.data., sep="-", into=c("cells","num"),
                             remove=F)
cells_ctl <- cells2 %>% dplyr::filter(scdata.meta.data.sampletype == "CTL") %>% 
  dplyr::select(rownames.scdata.meta.data., scdata.meta.data.predicted.id)
cells_NP <- cells2 %>% dplyr::filter(scdata.meta.data.sampletype == "NP") %>% 
  mutate(cell0 = paste0(cells, "-1")) %>% 
  dplyr::select(cell0, scdata.meta.data.predicted.id)
write.table(cells_ctl, "out/cells_ctl.txt", row.names = F, col.names = F,
            sep = "\t", quote = F)
write.table(cells_NP, "out/cells_NP.txt", row.names = F, col.names = F,
            sep = "\t", quote = F)
write.table(cells, "out/cells.txt", row.names = F, col.names = F,
            sep = "\t", quote = F)
cond <- data.frame(rownames(scdata@meta.data),scdata@meta.data$sampletype,
                   scdata@meta.data$predicted.id)
cond <- cond %>% mutate(cond = paste0(scdata.meta.data.sampletype, "_",
                                      scdata.meta.data.predicted.id )) %>%
  dplyr::select(rownames.scdata.meta.data.,cond)
write.table(con, "out/cond.txt", row.names = F, col.names = F,
            sep = "\t", quote = F)
con <- read.table("out/cond.txt")
head(con)
con$V2 <- factor(con$V2, levels = unique(con$V2),
                 labels = c("CTL_Cluster2", "NP_Cluster4", "NP_Cluster2",
                            "NP_Cluster5", "CTL_Cluster1",  "CTL_Cluster4",
                            "NP_Cluster1", "CTL_Cluster3",  "NP_Cluster3",
                            "NP_Cluster5"))
## bigwigs split
library(trackViewer)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggpubr)


f <- "Il10"
bw <- list.files(path= "data/aggr_multiome/subset_gex",
                 pattern= "Cluster*")

cols2 <- c("#00B0F6", "#00BF7D", "#F8766D", "#A3A500", "#E76BF3")
DefaultAssay(scdata) <- "ATAC"
Idents(scdata) <- "predicted.id"

bwig <- function(bw,w,tec,dir){
  pl <-  BigwigTrack(
    region = agr[[1]],
    bigwig = paste0("data/aggr_multiome/", dir, bw[w]),
    smooth = 3000,
    type = "coverage",
    y_label = tec,
    max.downsample = 3000, 
    downsample.rate = 0.1) +
    ggtitle(strsplit(bw[w], ".", fixed=T)[1][[1]][1]) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title.y = element_text(size = 9)) +
    xlab("")
  pl$layers[[1]]$aes_params$fill <- cols[w]
  return(pl)
}
pdf("figs/coverage_tracks.pdf", width = 16, height = 25)
for (f in features[-22]) {
  id <- get(f, org.SYMBOL2EG,)
  gr <- genes(TxDb, single.strand.genes.only=FALSE)[id]
  agr <- if (width(gr)[[1]]<120000) {
    resize(gr, 100000, fix = "center")
  } else {
    resize(gr, width(gr)[[1]]+100000, fix = "center")
  }
  
  # genesinrange <- mapRangesToIds(TxDb, agr, type = "gene", ignore.strand = T)
  # tracks <- sapply(genesinrange[[1]][[1]],
  #                  function(z) {
  #                    track <- geneTrack(z,TxDb)[[1]]
  #                    return(track)
  #                  })
  
  
  rn <- list()
  rn_ctl <- list()
  rn_t <- list()
  ac <- list()
  ac_ctl <- list()
  ac_t <- list()
  for (w in 1:length(bw)) {
    rn[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/")
    ac[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/")
    rn_ctl[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/control/" )
    ac_ctl[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/control/")
    rn_t[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/treated/")
    ac_t[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/treated/")
  }
  #library(ggpubr)
  ann <- AnnotationPlot(scdata,region=GRangesToString(agr))
  ann$layers[[4]]$aes_params$size <- 6
  plt <- ggarrange(rn[[1]], ac[[1]], rn[[2]], ac[[2]],
                   rn[[3]], ac[[3]], rn[[4]], ac[[4]],
                   rn[[5]], ac[[5]],
                   ann, ncol=1)
  plt <- annotate_figure(plt,
                         top = text_grob("CD1d and CTL", color = "grey",
                                         face = "bold", size = 16))
  plt_ctl <- ggarrange(rn_ctl[[1]], ac_ctl[[1]], rn_ctl[[2]], ac_ctl[[2]],
                       rn_ctl[[3]], ac_ctl[[3]], rn_ctl[[4]], ac_ctl[[4]],
                       rn_ctl[[5]], ac_ctl[[5]], ann, ncol=1)
  plt_ctl <- annotate_figure(plt_ctl,
                             top = text_grob("CTL", color = "grey",
                                             face = "bold", size = 16))
  plt_t <- ggarrange(rn_t[[1]], ac_t[[1]], rn_t[[2]], ac_t[[2]],
                     rn_t[[3]], ac_t[[3]], rn_t[[4]], ac_t[[4]],
                     rn_t[[5]], ac_t[[5]], ann, ncol=1)
  plt_t <- annotate_figure(plt_t,
                           top = text_grob("CD1d", color = "grey",
                                           face = "bold", size = 16))
  plt_all <- ggarrange(plt_ctl, plt_t, plt, ncol = 3)
  plt_all <- annotate_figure(plt_all,
                             top = text_grob(f, color = "black",
                                             face = "bold", size = 19)) 
  print(plt_all)
}
dev.off()


ac_5 <- markers_atac %>% filter(cluster == "Cluster5") %>%
  filter(avg_log2FC < 0)
ac_5gene <- closest$gene_name[closest$query_region %in% ac_5$gene] %>% unique()
rna_5_down <- markers_rna %>% filter(cluster == "Cluster5") %>%
  filter(avg_log2FC < 0)
rna_5_up <- markers_rna %>% filter(cluster == "Cluster5") %>%
  filter(avg_log2FC > 0)

writeClipboard(ac_5gene)
          
Vol <- data.frame(log2FC= rna_5$avg_log2FC,
                  sig= -log(rna_5$p_val_adj, 20), 
                  row.names = rna_5$gene)
Vol$S <- 0
Vol[which(Vol$sig<1.5), "S"] <- "NS"
Vol[which(Vol$log2FC>0.5 & Vol$S!="NS"), "S"] <- "UP"
Vol[which(Vol$log2FC<(-0.5) & Vol$S!="NS"), "S"] <- "DOWN"
Vol[which(Vol$S==0), "S"] <- "NS"
Vol$S <- as.factor(Vol$S)
-log(1.5,20)

library(ggrepel)
volcano <- ggplot(Vol, aes(x=log2FC, y=sig, color = S)) +
  geom_point (size = 1.4,
              show.legend = F) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 1, size = 0.3, col = "grey20") +
  geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
  scale_color_manual(values = c("skyblue2", "grey60", "green3")) +
  theme_light() +
  titol("Cluster 5 DE genes out of Diff OCR") +
  xlab("Fold Change") + 
  ylab("-log20 pval_adj") +
  scale_y_continuous(breaks = seq(0, , 50),
                     limits = c(-20, 400), expand = expansion(mult=0, add= c(1,0))) +
  scale_x_continuous(breaks = c(seq(-30,30, by = 1)),
                     labels = c(2^abs(-30:30)*sign(-30:30)),
                     guide = guide_axis(check.overlap = T)) +
  theme(axis.line = element_line(size = 0.3, colour = "grey20", linetype=1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  geom_text_repel(data= subset(Vol, Vol$S != "NS"),
                  aes(label = rownames(subset(Vol, Vol$S != "NS"))),
                  color = "black",
                  size = 4,
                  vjust = -0.1, max.overlaps = 20)
volcano
