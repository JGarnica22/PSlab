# Analysis of V(D)J, Gene Expression and Feature Barcode single-cell data :art:
The 5' Chromium Single Cell Immune Profiling Solution with Feature Barcode technology enables simultaneous profiling of V(D)J repertoire, cell surface protein, antigen specificity and gene expression data. This folder contain the tool and indications necessary to analyze together all this data.
</br>

As in other single-cell pipelines guides, there is a first part run in shell based on cellranger pipelines which basically involves alignment, filtering, barcode counting and TCR or BCR aligment. These should preferably be run in the cluster and a loop to run all cellranger steps for all the samples can be found in [cellranger_multi_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/cellranger_multi_loop.sh). After this, this guideline recommends downloading cellranger outputs and importing them into Seurat package to do the downstream analysis, visualization and TCR clonotype integration. A script with Seurat functions to use is found at [10X_RNA_ADT_VDJ.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/10X_RNA_ADT_VDJ.R).

Go through all these steps with [10X_multi_VDJ_RNA_ADT.md](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/10X_multi_VDJ_RNA_ADT.md). For more information about cellranger-arc also visit [10X webpage](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi) or [Satija lab](https://satijalab.org/). 
</br>
