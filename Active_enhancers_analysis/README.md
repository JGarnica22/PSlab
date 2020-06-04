# Active enhancers analysis :signal_strength:

These scripts can be used to infer the localization of active enhancers based on **ChIP** and **ATAC-seq data**. Active regulatory regions can be obtained by overlapping H3K27ac-ChIP and ATACseq data, as active regulatory regions are expected to be both open and histone-acetylated. To obtain active enhancers, promoters are excluded since they share these characteristics with active enhancers. 

You will need the following input data files:

*NOTE: Bear in mind that scripts are prepared to look for these files in your directory (/data), and they should include the name of the technique (ATAC, ChIP...) and the name of the population or sample as filled in the script.*
- ChIP annotated peaks for each of your samples
- ATACseq annotated peaks for each of your samples
- ATACseq shared regions (OCRs) between samples
- DMR data file (Differentially Methylated regions) comparing your samples
- DESeq2 comparing your populations, as a data frame detailed in [DE_analysis_RNAseq](https://github.com/patriciasolesanchez/PSlab/blob/master/DE_analysis_RNAseq/DE_analysis_RNAseq_1vs1.R).

In this folder you will find the following scripts:

* [Active_enhancers_all.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/Active_enhancers_all.R) is a comprehensive script overlapping ChIP, ATAC, DMR and DESeq2 data to compare active enhancers between samples, providing outputs files for each analysis and a final **summary table** and **bar chart** quantifying the differences between populations.

After active enhancers localization, the script predicts possible genes that may be targets of these active enhancers: it annotates each enhancer to genes found 100 Kb around these active enhancers. This step in this script needs to be performed using **Terminal** `beedtools window` tool taking few seconds, as it is detailed in the same script. 

* Alternatively, if you can <ins>only work with R</ins> (can't use Terminal), see [`Active_enhancers_all_R.R`](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/Active_enhancers_all_R.R). This script performs this same annotation step (search for genes around active enhancers) without using Terminal, as a step integrated in the whole R loop. However, **it takes several minutes to complete**.

* Finally, [`Genomic_screening_around_gene.R`](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/Genomic_screening_around_gene.R) is useful to associate data to a particular gene (or genes). It allows you to look for certain genomic elements (e.g. ATAC OCRs, ChIP peaks, active enhancers, etc.) around a certain gene. Basically, it creates a loop iterating all the genes of interest and looks for elements around these genes by overlapping with `GRanges` objects. Finally, it returns a data frame for all the localizations of elements regions close to the gene.
