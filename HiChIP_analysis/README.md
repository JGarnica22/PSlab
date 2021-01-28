# HiChIP analysis :thread:

This pipeline can be used to analyse HiChIP data. Briefly, HiChIP is a technique used to identify regions that are in contact in the 3D structure of the genome. In particular, we are connecting promoters to the regions they are contacting with by using a promoter enrichment procedure in the HiChIP protocol (see reference in /doc).

HiChIP data can be analysed using the 'HiC-Pro' pipeline from Servant (for more information see the [documentation](http://nservant.github.io/HiC-Pro/QUICKSTART.html) or their [GitHub page](https://github.com/nservant/HiC-Pro)).

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. HiC-Pro supports the main Hi-C protocols, including digestion protocols as well as protocols that do not require restriction enzyme such as DNase Hi-C. In practice, HiC-Pro can be used to process dilution Hi-C, in situ Hi-C, DNase Hi-C, Micro-C, capture-C, capture Hi-C or HiChIP data. The pipeline is flexible, scalable and optimized. HiC-Pro is sequential and each step of the workflow can be run independantly.

HiC-Pro requires several dependencies; we are using Singurality to install the HiC-Pro pipeline and all its dependencies at once.  
</br>

In order to run HiC-Pro, you need to:

1) Prepare annotation files.

2) Put all your input files in a rawdata folder.

3) Set up your configuration file _config-install.txt_.

4) Run HiC-Pro.  
</br>

## Annotation files


## Raw data folder


## Configuration file


## Run HiC-Pro

´´´´
MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE
´´´´


You can run your HiC-Pro analysis in the cluster by directly submitting the _hicpro.sh_ script.





