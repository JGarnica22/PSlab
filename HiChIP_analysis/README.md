# HiChIP analysis :thread:

This pipeline should be used to analyse HiChIP data. Briefly, HiChIP is a technique used to identify regions that are in contact in the 3D structure of the genome. In particular, we are connecting promoters to the regions they are connected to by using a promoter enrichment procedure in the HiChIP protocol (see reference in /doc).

HiChIP data can be analysed using the 'HiC-Pro' pipeline from Servant (for more information see [documentation](http://nservant.github.io/HiC-Pro/QUICKSTART.html) or their [GitHub page](https://github.com/nservant/HiC-Pro)).

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. HiC-Pro supports the main Hi-C protocols, including digestion protocols as well as protocols that do not require restriction enzyme such as DNase Hi-C. In practice, HiC-Pro can be used to process dilution Hi-C, in situ Hi-C, DNase Hi-C, Micro-C, capture-C, capture Hi-C or HiChIP data. The pipeline is flexible, scalable and optimized. HiC-Pro is sequential and each step of the workflow can be run independantly.

In order to run HiC-Pro, which requires several dependencies, we are using Singurality.

Before running HiC-Pro, you need to:

1) Set up your configuration file 
