# Single cell RNAseq analysis :zap:

Chromium Single Cell data (10x data) can be analyzed using Cell Ranger and Seurat. First, the analysis pipeline in Cell Ranger performs sample demultiplexing, barcode processing, and single cell 3' gene counting. Then, gene count matrices can be used in Seurat to perform clustering and differential expression analysis.

## Cell Ranger 
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
[Josep]


## Seurat

Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data. Seurat is developed and maintained by the [Satija lab] (https://satijalab.org/seurat/). For guided tutorials on how to use Seurat, visit this [website] (https://satijalab.org/seurat/vignettes.html).

The input for Seurat is the output of Cell Ranger, basically an aggr ("aggregate") file that contains 3 files:  
     barcodes.tsv.gz  
     features.tsv.gz  
     matrix.mtx.gz  


