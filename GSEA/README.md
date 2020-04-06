# GSEA ANALYSIS :octocat:
Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states (e.g. phenotypes).

The R package used, `fgsea`, implements an algorithm for fast gene set enrichment analysis. Using the fast algorithm allows to make more permutations and get more fine grained p-values, which allows to use accurate stantard approaches to multiple hypothesis correction.

Analysis is performed with differential expression data from previous analysis, e.g. from DESeq2 data.

For this analysis you will need:
  1. Differential expression data file containing genes and log2FoldChange
  2. Gene sets of interest, i.e. gmt file from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). We provide in this GitHub repository with the Hallmark gene sets (Classification of coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes). If you are interested in studying enrichment for other gene sets, go to website and dowloand your gmt
**IMPORTANT! It needs to contain genes in sets labelled as Gene symbols, do not use Entrez IDs.** 
**Gene sets from MSigDB are made with human genes, if using mouse data, transform gmt.file to mouse gene names.**
   
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your differential expression data file in the read.table
  4. Indicate populations (conditions) in the comparison
  5. Indicate gmt.file (Gene sets) to use for the analysis
  
Run all the script :smile:
