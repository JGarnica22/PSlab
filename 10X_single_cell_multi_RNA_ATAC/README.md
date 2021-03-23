# Single cell Multiome ATAC + Gene Expression analysis :dizzy:
In this folder you can find the tools and steps needed to analyse Chromium Single Cell Multiome ATAC + Gene Expression sequencing data. This will allow us to study both gene expression and chromatin accessibility, and their relationship, since both measurements are taken on the very same cell.

As in other single-cell pipelines guides, there is a first part run in shell based on cellranger pipelines which basically involves alignment, filtering, barcode counting and peak calling. These should preferably be run in the cluster and a loop to run all cellranger steps for all the samples can be found in [cellranger_arc_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/cellranger_arc_loop.sh).
After this, this guideline recommends downloading cellranger outputs and importing them into Seurat package to do the downstream analysis and visualization. A script with Seurat functions to use is found at [10X_ATAC_RNA.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/10X_ATAC_RNA.R).

Go through all these steps with [single_cell_multi_ATAC_RNA.md](https://github.com/patriciasolesanchez/PSlab/blob/master/10X_single_cell_multi_RNA_ATAC/single_cell_multi_ATAC_RNA.md). For more information about cellranger-arc also visit [10X webpage](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc) or [Satija lab](https://satijalab.org/). 
</br>
