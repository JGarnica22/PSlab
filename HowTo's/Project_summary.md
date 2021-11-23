# Project data summary
* Josep Garnica

* Date of last modification: 23/11/21


## Overview
This analysis is based on ChIP-seq data obtained from SANTAMARIA_16 project on September 2020. It included 3 samples for K27Ac, 3 for K4me3 and 2 inputs samples. More information at [sequencing projects](https://docs.google.com/spreadsheets/d/1lBpYgDUl5U1M_560ZxvkfJOyA9f4P76L/edit#gid=1548358353) on drive.

## Data
Data for this analysis involve bam files and peak datasets (.xlsx). These were obtained using bowtie2 and macs2, as described in the file `bowtie2_macs.sh` file (*download and save in the folder the bash script too*!).

## Docs
* `SANTAMARIA_16.xlsx`: sample list and association with `.fastq` files.

## Out
Output files:
* Peak dataset of each sample.
* Differential peakdata set for each comparison, indicated in the file tittle.
* Shared peaks between samples.

## Figs
Output figures, comparisons specicifies in file tittles:
* Barplot of peak datasets
* Cross-corellation heatmap
* Volcano plot

## Scripts
* Bash scripts for aligment and peak calling: `bowties_macs.sh`
* R script for differential binding analysis: `DiffBind.R`
* R script for active enhancers analysis: `Active_enhancers.R`
    * Overlap of differential ATAC peaks with K27Ac peaks.
* R script for active enhancers using another packages: `Active_enhancers_2.R`
