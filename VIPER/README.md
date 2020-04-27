# VIPER ANALYSIS :octocat:
This script can be used to run **Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER)** on differential expression analysis data. The `VIPER` algorithm allows computational inference of protein activity from gene expression data. It uses the expression of genes that are most directly regulated by a given protein, such as the targets of a transcription factor (TF), as an accurate reporter of its activity.

VIPER is run on a previous network defining interactions. In our case, we will infer TF activity from a TF-target network that was generated with **ARACNe (Algorithm for the Reconstruction of Accurate Cellular Networks)**. 

**IMPORTANT!! In this script we are providing a specific cellular network created on mouse CD4 T cell data.** If you want to use `VIPER` on a different cellular network you should either download it from an available source or create your own with `ARACNe`.<br/>
<br/>

For the script you will need the following input files:
  1. Differential expression data file containing genes and log2FoldChange, for example, from `DESeq2` data.
  2. rlog.norm.counts data file from DESeq2 script.
  3. Normalized counts matrix used to build your network (tissue-specific gene expression data that was used to run `ARACNe`). Here we provide a normalized counts matrix for CD4 T cell samples (data was downloaded from GEO). Find _norm.matrix.CD4.txt_ in data folder.
  4. Cellular network defining interactins (ARACNe output). Here we provide a mouse CD4 T cell network defining targets for a list of TFs. Find _mouse_CD4_T_cell_network_ARACNe.txt_ in data folder.<br/>
 <br/>
   
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your differential expression data file in the read.table
  4. Indicate populations (conditions) in the comparison
  5. Run all the script :smiley: <br/>
 <br/>

You will generate the following outputs:<br/>
1. Histogram with regulon distribution
2. Table of TF activity
3. Volcano plot
