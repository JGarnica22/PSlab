# DoRothEA :octocat:

This script can be used to infer **transcription factor (TF) activity** based on target differential expression. It is run using the `DoRothEA` regulon and the Virtual Inference of Protein-activity by Enriched Regulon analysis (`VIPER`) tool. <br/>

`DoRothEA` is a consensus regulon that links thousands of TFs to their targets based on different types of evidences: (i) literature curated resources, (ii) interactions derived from ChIP-seq data, (iii) in silico prediction of 
TF binding on motifs and (iv) reverse-engineered regulons from large gene expression datasets. <br/> 

Each interaction between TF and target gene is given a confidence from A to E. Interactions that are supported by all the four sources of evidence, manually curated by experts in specific reviews or supported both in at least two curated resource are considered to be highly reliable and were assigned an A score. Scores B-D are reserved for curated and/or ChIP-seq interactions with different levels of additional evidence. Finally, E score is used for interactions that are uniquely supported by computational predictions.<br/>

You can find more info in: https://saezlab.github.io/dorothea/ <br/>
<br />

For the script you will need the following input files: <br/>
  1. rlog normalized counts and differential expression data file from RNAseq analysis. <br/>
  2. Regulon for the species of the sample (human or mouse) as _.rdata_. You can download it from this repository in data folder. <br/>
  <br />

To do the analysis follow these steps: <br/>
  1. Open R script
  2. Change working directory
  3. Indicate the exact name of your data files <br/>
  4. Indicate populations (conditions) in the comparison <br/>
  5. Indicate which species you are working with: "human" or "mouse"
  6. Indicate confidence you want to use to do your analysis (only A, or only A,B,C, etc).<br/>
  7. Run all the script :smiley:<br />
  <br />

You will generate the following outputs:<br/>
1. TF activity inferences based on normalized gene expression of all confidences <br/>
2. TF activity inferences based on normalized gene expression of selected confidences <br/>
3. TF activity inferences based on gene expression signature (GES) of all confidences <br/>
4. TF activity inferences based on gene expression signature (GES) of all confidences with high significances <br/>

and the following figures: <br/>
1. Heatmap of TF activity inferences based on normalized gene expression of all confidences <br/>
2. Heatmap of TF activity inferences based on normalized gene expression of selected confidences <br/>
3. Volcano plot of TF activities differentially infered <br/>
4. Graphs plot with significantly differentially expressed TF in each confidence group <br/>
