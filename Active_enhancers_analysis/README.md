# Active enhancers analysis :signal_strength:

This script allow to infer the localization of active enhancers based on **ChIP** and **ATAC-seq data**. Active enhancer localization is inferred by overlapping ChIP and ATACseq data, as active enhancer chromatin regions are expected to be both open and histone-acetylated. Then promoters are excluded since they share these characteristics with active enhancers. 

You will need the following input data files:

*NOTE: Bear in mind that scripts are prepared to look for these files in your directory (/data), so these files should also include the name of the technique (ATAC, ChIP...) and the name of the population or sample as filled in the script.*
- ChIP annotated peaks for each of your samples
- ATACseq annotated peaks for each of your samples
- ATACseq shared OCR between samples
- DMR data file (Differentially Methylated egions) comparing your samples
- DESeq2 comparing your populations, as a data frame detailed in [DE_analysis_RNAseq](https://github.com/patriciasolesanchez/PSlab/blob/master/DE_analysis_RNAseq/DE_analysis_RNAseq_1vs1.R).

[`Active_enhancers_all.R`](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/Active_enhancers_all.R) is a comprehensive script overlapping ChIP, ATAC, DMR and DESeq2 data to compare active enhancers between samples, giving as outputs files for each analysis and a final **summary table** and **bar chart** for all of them.

After active enhancers localization, the script predict possible genes that may affected by these active enhancers just by proximity, so it returns the genes 100 Kb around these active enhancers. This step in this script need to be performed using **Terminal** `beedtools window` tool taking few seconds, as it is detailed in the same script. 

Alternatively, [`Active_enhancers_all_R.R`](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/active_enhancers_all_R.R) script allow the search of these genes around active enhancers without using Terminal which is integrated in the whole loop, however, it takes several minutes to complete.

Finally, [`Genomic_screening_around_gene.R`](https://github.com/patriciasolesanchez/PSlab/blob/master/Active_enhancers_analysis/Genomic_screening_around_gene.R) allow to look for certain genomic elements, such as active enhancers but any other regions of interest, around a certain gene or other element. Basically, it creates a loop iterating all the genes of interest and looks for elements around these genes by overlapping with `GRanges` objects. Finally, it returns a data frame for all the localizations of active enhancers or other regions close to the gene.
