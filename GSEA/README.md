# GSEA ANALYSIS :octocat:
Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically
significant, concordant differences between two biological states
(e.g. phenotypes).

The R package used, fgsea, implements an algorithm for fast gene set enrichment analysis. 
Using the fast algorithm allows to make more permutations and get more fine grained p-values, which allows to use accurate stantard approaches to multiple hypothesis correction.

Analysis is performed with differential expression data from previous analysis, for example from DESeq2 data.

For this analysis you will need the following file:
  1. Differential expression data file containing genes and log2FoldChange.
   
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your differential expression data file in the read.table
  
 
 Run all the script
