# Dorothea

This script infers trasncription factor (TF) activity based on target differential expression. This based on different sources of data <br/>
including (i) manually curated repositories, (ii) interactions derived from ChIP-seq binding data, (iii) in silico prediction of <br/>
TF binding on gene promoters and (iv) reverse-engineered regulons from large gene expression datasets. <br/>
Each inference between TF and target gene is given a confidence from A to E. Interactions that are supported by all the four sources <br/>
of evidence, manually curated by experts in specific reviews or supported both in at least two curated resource are considered to <br/>
be highly reliable and were assigned an A score. Scores B-D are reserved for curated and/or ChIP-seq interactions with different <br/>
levels of additional evidence. Finally, E score is used for interactions that are uniquely supported by computational predictions.<br/>
<br/>
For this script you will need: <br/>
  1. rlog normalized counts and differntial expression from RNAseq analysis, or at least the raw counts. <br/>
  2. Regulon for the species of the sample (human or mouse) as Rdata or csv. <br/>
<br/>
As you run the script you will need to indicate:<br/>
  1. The name of the required files. When indicating DE_matrix, you need to choose which comparison in case more than one sample or
  condition is present. <br/>
  2. Populations or conditions of your analysis. <br/>
  3. Confidence you want to use to do your analysis (only A, or only A,B,C, etc).<br/>
<br/>
