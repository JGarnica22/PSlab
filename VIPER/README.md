# VIPER ANALYSIS :octocat:
This script can be used to perform **VIPER** on differential expression analysis data. 
The VIPER (Virtual Inference of Protein-activity by Enriched Regulon analysis) algorithm[12] allows computational inference of protein activity, on an individual sample basis, from gene expression proﬁle data. It uses the expression of genes that are most directly regulated by a given protein, such as the targets of a transcription factor (TF), as an accurate reporter of its activity. 

The R package used, `viper`, implements an algorithm that inference protein activity from gene expression data

For the script you will need the following input file:
  1. Differential expression data file containing genes and log2FoldChange, for example, from `DESeq2` data.
  2. rlog.norm.counts data file from DESeq2 script.
  3. Txt file containing normalized matrix for your type of samples (tissue-specific data needs to be downloaded from GEO).
  4. Txt file containing ARACNe output.
   
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your differential expression data file in the read.table
  4. Indicate populations (conditions) in the comparison
   
Run all the script :smiley:

You will generate the following outputs:<br/>
1. Histogram with regulon distribution
2. Table of TF activity
3. Vocano plot

