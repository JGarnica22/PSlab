# DE_analysis_RNAseq

These scripts can be used to perform differential expression analysis (DE) on RNAseq data using DESeq2 package from Bioconductor. The package estimates variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.<br />

You will find one script to run DE analysis on 2 conditions, _DE_analysis_1vs1.R_. If you have several conditions, you can use this script and make all comparisons you want independently, running them 1 at a time, or you can use the script _DE_analysis_multiple_samples.R_ to compare in one run all conditions to the control. Using this method you will get less significant genes, because p-value correction on multiple comparisons is more strict. We recommend you compare populations 1 vs 1.<br />


For the script(s) you will need the following input file:
1. Raw counts data file.

To do the analysis follow these steps:
1. Open R script
2. Change working directory
3. Indicate the exact name of your raw counts file as "expfile"
4. Indicate the population(s) to compare. *_Note that always control population needs to be indicated first_
5. Indicate which species you are working with: "human" or "mouse"

Run all the script 😄

You will generate the following outputs:
1. rlog-normalized counts
2. PCA (Principal component analysis)
3. DESeq2 DE analysis
4. Volcano plot
5. Heatmap
