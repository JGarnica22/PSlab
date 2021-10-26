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
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Huinkt/")

# Read the data:
# The Read10X function reads in the output of the cellranger pipeline from 10X 
# (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
# returning a unique molecular identified (UMI) count matrix. The values in this 
# matrix represent the number of molecules for each feature (i.e. gene; row) that 
# are detected in each cell (column).

# Indicate in the Read10X function, the directory of your input files (3 files must be
# found in the data.dir: matrix.mt.gz, barcodes.tsv.gz and features.tsv.gz)
hu_data <- Read10X(data.dir = paste0(getwd(),"/data/human"))
mo_data <- Read10X(data.dir = paste0(getwd(),"/data/mouse"))

#get orthologs genes
getBM(attributes = c("mmusculus_homolog_associated_gene_name", "mmusculus_homolog_chromosome"),
      filters = "hgnc_symbol",
      values = c("E2F3"),
      mart = human.ensembl)


mo_ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
hu_ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

bmmix <- getLDS(attributes = "hgnc_symbol",
                filters = "hgnc_symbol", values = row.names(hu_data), mart = hu_ensembl,
                attributesL = "mgi_symbol" , martL = mo_ensembl)


## Cross-species data integration, change human gens to mouse genes, human genes without an homologous gene will be deleted!

z  <- as.data.frame(rownames(hu_data))
names(z) <- "HGNC.symbol"
z$id <- 1:nrow(z)
y <- merge(z, bmmix, by = "HGNC.symbol" , all = T)
y <- y[order(y$id),]
hu_data@Dimnames[[3]] <- y$MGI.symbol
hom_data <- subset(hu_data, is.na(hu_data@Dimnames[[3]]) == F)

# We next use the count matrix to create a Seurat object. The object serves as a 
# container that contains both data (like the count matrix) and analysis (like PCA, 
# or clustering results) for a single-cell dataset.
hu_scdata <- CreateSeuratObject(counts = hu_data, min.cells = 3, min.features = 200)
mo_scdata <- CreateSeuratObject(counts = mo_data, min.cells = 3, min.features = 200)
