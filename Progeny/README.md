# PROGENY ANALYSIS :octocat:
Pathway RespOnsive GENes (PROGENy) is a bioinformatic tool to infer pathway activity based on transcriptomic data.  
Based on perturbation experiments, it asigns scores to each gene involved in a pathway depending on how correlates each gene to each pathway activation. Then it calculates pathway activity based on expression level/deregulation of these genes in your sample(s).

For more information go to: https://saezlab.github.io/progeny/ 

Analysis can be performed both on raw counts data to assess pathway activity in single sample (followed by a linear regression model to study difference between conditions) or on differential expression data (e.g. from DESeq2).

For this analysis you will need beforehand the following files:
  1. Raw counts data file.
  2. Progeny matrix for the species you are working with (human or mouse), available to download in this repository, version 2020.
  3. DESeq2 analysis data file (it is possible to do it inline here too).
  
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your raw counts file as "expfile"
  4. Indicate the exact name of your DESeq2 analysis data file as "Desq2file"
  5. Indicate which species you are working with: "human" or "mouse"  
 
Run all the script :smiling:
