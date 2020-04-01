# PROGENY ANALYSIS :octocat:
Pathway RespOnsive GENes (PROGENy) is a bioinformatic tool to infer pathway activity based on transcriptomic data.  
Based on experiments, it gives scores to each gene involved in a pathway depending on how correlates each gene to each pathway activation and then infer pathway activity based on upregulation of these genes in your sample(s).

Analysis is performed both on raw counts data to assess how behave each replicate and also on DESeq2 data.

For this analysis you will need beforehand the following files:
  1. Raw counts data file.
  2. Progeny matrix for the species you are working with (human or mouse), available to download in these repository, version 2020.
  3. Desq2 analysis data file (it is possible to do it inline here too).
  
To do the analysis follow these steps:
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your raw counts file as "expfile"
  4. Indicate the exact name of your DESeq2 analysis data file as "Desq2file"
  5. Indicate with which species you are working: "human" or "mouse"  
 
 Run all the script