# PROGENY ANALYSIS :octocat:
This script can be used to perform **Pathway RespOnsive GENes (PROGENy) analysis**. `PROGENy` is a bioinformatic tool to infer pathway activity based on transcriptomic data. Based on perturbation experiments, it asigns scores to each gene involved in a pathway depending on how correlates each gene to each pathway activation. Then it calculates pathway activity based on expression level/deregulation of these genes in your sample(s).

For more information go to: https://saezlab.github.io/progeny/ 

Analysis can be performed both on raw counts data to assess pathway activity in single sample (followed by a linear regression model to study difference between conditions) or on differential expression data (e.g. from `DESeq2`).

For the script you will need the following input file(s):
  1. Raw counts data file.
  2. Progeny matrix for the species you are working with (human or mouse), available to download in this repository, version 2020.
  3. DESeq2 analysis data file (it is possible to do it inline here too).
  
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your raw counts file as "expfile"
  4. Indicate the exact name of your DE analysis data file as "Desq2file"
  5. Indicate which species you are working with: "human" or "mouse"  
 
Run all the script :smile:

You will generate the following outputs:<br/>
1. Linear regression of progeny permutations between samples or conditions<br/>
2. Progeny scores based on DESeq2 results<br/>

and the following figures:<br/>
1. Heatmap with the progeny permutations between samples or conditions<br/>
2. Dotplot graph of differences in score between samples for each pathway<br/>
3. Graphbars with contrast score for each progeny pathway based on DESeq2 results<br/>
4. Dotplot with genes most differentialy expressed within each pathway<br/>

