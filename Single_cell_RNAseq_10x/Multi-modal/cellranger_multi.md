# Analyzing V(D)J, Gene Expression and Feature Barcode with cellranger multi :octopus:


The 5' Chromium Single Cell Immune Profiling Solution with Feature Barcode technology enables simultaneous profiling of V(D)J repertoire, cell surface protein, antigen specificity and gene expression data. The cellranger multi pipeline enables the analysis of these multiple library types together. The advantage of using the multi pipeline (as opposed to using `cellranger vdj` and `cellranger count` separately) is that it enables more consistent cell calling between the V(D)J and gene expression data.

**Note**: These guide lines assume that you already have fastq files either recieved from the facility or produced with `cellranger mkfastq`.

For more information visit [10X genomics webpage](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi).
