# This script shows how to load data from 10X multi-modal experiments including
# gene expression, antibody capture (or ADT: Antibody-derived tags) and VDJ data.

# In this script it is assumed previous knowledge on Seurat package and non-multi-modal 
# data analysis. Here only the differences in multi-modal analysis will be shown. 
# For more info visit 
# https://satijalab.org/seurat/articles/multimodal_vignette.html

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)


# Seurat is able to analyze data from multimodal 10X experiments, in this case
# gene expression and ADT data are comprised in the same file, while VDJ data
# is found in another one and must be added as metada to the Seurat object.

# Load 10X data as usual, you will be informed that this is a list of 2
# including gene expression and ADT data.
data <- Read10X(data.dir = "data/count/filtered_feature_bc_matrix")

# Firstly, create a Seurat object with gene expression data
scdata <- CreateSeuratObject(counts = data[["Gene Expression"]], 
                             min.cells = 3, min.features = 200)

# Subset your data and add the metadata (samples, condition...) and
# normalize the data if needed
scdata <- NormalizeData(scdata, normalization.method = "LogNormalize",
                        scale.factor = 10000) #only if not normalized in cellranger!!


# After that you can add the ADT data to the Seurat object
scdata[["ADT"]] <-  CreateAssayObject(data[["Antibody Capture"]][,colnames(scdata),
                                      drop = F])

# Again normalize the ADT data if not previously
scdata <- NormalizeData(scdata, assay = "TET", normalization.method = "CLR")


# Add VDJ data, in this example TCR.
# VDJ data cannot be included as an assay and must be added as metadata, as
# is categorical data. Then clonotypes can be filtered and represented accordingly.

# Import vdj data coming from cellranger pipeline, saved in vdj_t directory
tcr <- read.csv("data/vdj_t/filtered_contig_annotations.csv", sep=",")
# This file gives us the productive alpha and beta chain of each cell
# In order to analyze the data we need each cell to appear only once,
# so we will be using the clonotypes:

# Remove the -1 at the end of each barcode.
# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have the same clonotype.
tcr$barcode <- gsub("-1", "", tcr$barcode)
tcr <- tcr[!duplicated(tcr$barcode), ]

# Only keep the barcode and clonotype columns
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode", "raw_clonotype_id")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

# Load the clonotype-centric info
clono <- read.csv("data/vdj_t/clonotypes.csv", sep=",")

# Slap the nt sequences onto our original table by clonotype_id
# you can add more columns if you want this information
tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_nt")])

# Set barcodes as rownames
rownames(tcr) <- tcr[,"barcode"]
tcr[,"barcode"] <- NULL

# Add to the Seurat object's metadata
scdata <- AddMetaData(scdata, metadata=tcr)


# At this point you have your Seurat object with all the multi-modal
# data loaded and you can proceed with your analysis.

# Here we will show some functions unique and handy with this type of data

# Label cells based on clonotype id on tSNE plots:
DimPlot(scdata, group.by = "clonotype_id", reduction = "tsne", dims = c(1,2))
# Do the same as above but only show cells with productive clonotype assignation:
DimPlot(subset(scdata, subset = clonotype_id != "NA") , group.by = "clonotype_id",
        reduction = "tsne", dims = c(1,2))

# If you have multiple ADT variables you can do a scatter plot with those
# and even plot data from different assays

# Note that when working with different assays you should specify which
# you want to work with by:
DefaultAssay(scdata) <- "ADT"

# Alternately, we can use specific assay keys to specify a specific modality,
# thus we do not need to specify the default assay.
Key(cbmc[["RNA"]])
## [1] "rna_"
Key(cbmc[["ADT"]])
## [1] "adt_"

# Scatter for the same genes of different ??????
FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "rna_CD3E", pt.size = 1)

# Feature plot for an ADT
FeaturePlot(scdata, reduction = "tsne",
            features = "PE-TotalSeqC", 
            max.cutoff = 3,
            cols = c("grey", "red"), 
            combine = T)

