# GSEA ANALYSIS :octocat:
This script can be used to perform Gene Set Enrichment Analysis (GSEA) on differential expression analysis data. It is a computational method that determines whether an a priori defined set of genes shows statistically
significant, concordant differences between two biological states (e.g. phenotypes).

The R package used, fgsea, implements an algorithm for fast gene set enrichment analysis. Using the fast algorithm allows to make more permutations and get more fine grained p-values, which allows to use accurate stantard approaches to multiple hypothesis correction.

For the script you will need the following input file:
  1. Differential expression data file containing genes and log2FoldChange, for example, from DESeq2 data.
   
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your differential expression data file in the read.table
  4. Indicate which species you are working with: "human" or "mouse"
  
 
Run all the script :smiley:

You will generate the following outputs:
1. rlog-normalized counts
