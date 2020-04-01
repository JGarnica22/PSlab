#This script is to perform DE analysis with RNAseq data
#It was created for R 3.6.2 version (2019-03-30)
#Copyright (C) 2020  Patricia Sole Sanchez
# Check if required packages are installed, if not install:
cran.packages <- c("BiocManager","ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
                   "broom", "readr", "gem_point", "writexl")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}
bioc.packages <- c("DESeq2", "biomaRt", "progeny", "apeglm")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}
# Load packages:
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(writexl)
library(progeny)
library(magrittr) #To use command like %%
library(dplyr)
library(writexl)

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/Progeny")
# Specify parameters to be used along the script:
## Indicate name of txt file containing expression data (raw counts)
expfile <- "Partek_TFH_Raw_counts.txt"
##Indicate name of txt file containing DESeq2 analysis result previously performed
Desq2file <- "Deseq2_TFH_vs_Th0.txt"
## Indicate populations of interest, to be compared:
# First indicate control population, then sample:
pop <- c("Th0", "TFH")
#Indicate "human" or "mouse" experiment
species <- "mouse"
model <- read.table(paste0("data/progeny_matrix_",species,".txt"), header = TRUE, sep = "\t", dec = ".", quote = "") 
#Prepare the gene expression matrix for PROGENy:
counts <- read.table(paste0("data/",expfile),
                               header = T,
                               sep = "\t",
                               quote = "", 
                               dec = ".")
#Prepare the DESeq2 results
DESeq2 <- read.table (file = paste0("data/",Desq2file),
                      sep = "\t", 
                      quote = "",
                      dec = ".")

# Set gene_name as rows labels
row.names(counts) <- counts$gene_name
# Eliminate unnecessary columns, we only want gene_name and counts for the different populations.
counts2 <- counts[, c(grep(pop[1], colnames(counts)), grep(pop[2], colnames(counts)))]
# Change names of columns to shorter, easier ones:
names (counts2) <- c(paste0(rep(pop[1], length(grep(pop[1],colnames(counts2)))),
                             1:length(grep(pop[1],colnames(counts2)))), 
                      paste0(rep(pop[2], length(grep(pop[2],colnames(counts2)))),
                             1:length(grep(pop[2],colnames(counts2)))))
# Create conditions:
sample_info <- data.frame(cell_type = c(rep(pop[1], length(grep(pop[1], colnames(counts2)))), 
                                        rep(pop[2], length(grep(pop[2], colnames(counts2))))),
                          replicates = c(1:length(grep(pop[1],colnames(counts2))),
                                         1:length(grep(pop[2],colnames(counts2)))),
                          row.names = names (counts2))
# Transform the matrix by rounding the counts and transforming the values from "numeric" to "integer"
counts2_r <- apply(counts2, c(1,2), round)
counts2_i <- apply(counts2_r, c(1,2), as.integer)
# Generate the DESeq2 dataset:
dset = DESeqDataSetFromMatrix(countData = counts2_i,
                              colData= sample_info, 
                              design= ~ cell_type)
dset <- estimateSizeFactors(dset)
dset <- estimateDispersions(dset)
emat <- getVarianceStabilizedData(dset)

#We can not use the progeny function to score the expression matrix, because that was made for
#human. It will be done using a scoring matrix for mouse:
#Read the scoring matrix for progeny:
row.names(model) <- model[,1]
model <- model[,-c(1)]
#Select gene names that appear both in the expression matrix and in the scoring matrix:
common_genes <- intersect(rownames(emat), rownames(model))
#Filter both matrixes to only contain the genes that are shared
#1. Filter gene expression matrix:
emat_matched <- emat[common_genes, , drop = FALSE] %>%
  t()
#The option drop = FALSE ensures that the object will continue to be a data.frame
#Also the expression matrix should be transposed so that gene names are in rows and columns in each matrix
#That is a property of matrixes (to multiply, they should be one transposed from the other)
#2. Filter weight matrix:
model_matched <- model[common_genes, , drop = FALSE] %>%
  data.matrix()

#SANITY CHECK: gene names in the expression matrix should be on the columns and in the progeny
#scoring matrix they should be on the rows
stopifnot(names(emat_matched) == rownames(model_matched))

#Multiply the 2 matrixes to obtain the progeny scores:
progeny_scores = emat_matched %*% model_matched
#We will scale by pathways. For that, we use the function "scale", which directly scales by columns (pathways)
progeny_scores_scaled = scale(progeny_scores, scale = TRUE)

#Checking for differences between the groups:
#We check if the control is different to the treated condition using a linear model:
#Stablish control samples
controls <- sample_info$cell_type == pop[1]
#Perform linear regression to calculate enrichment of pathways:
result <- apply(t(progeny_scores_scaled), 1, function(x) {
  broom::tidy(lm(x ~ !controls)) %>%
    filter(term == "!controlsTRUE") %>%
    select(-term)
})
result <- mutate(Pathway=names(result), bind_rows(result))
results <- as.data.frame(result[,c(5,1:4)])
row.names(results)<-results$Pathway
results <- results[,-1]
write_xlsx(results, "output/Progeny_scores_scaled.xlsx")

#To check significance of progeny pathway scorings, we will perform Progeny on several (usually n=1000)
#permutations: use `runProgenyFast` function:
#'\code{runProgenyFast}
#'#'This function is designed to compute progeny pathway scores and assess their significance using a gene sampling based permutation strategy, 
#'for a series of experimental samples/contrasts.
#'#'@param df A data.frame of n*m+1 dimension, where n is the number of omic features to be considered and m is the number of samples/contrasts.
#'The first column should be the identifiers of the omic features. These identifiers must be coherent with the identifers of the weight matrix.
#'@param weight_matrix A progeny coeficient matrix. the first column should be the identifiers of the omic features, and should be coherent with the identifiers provided in df.
#'@param k The number of permutations to be preformed to generate the null-distribution used to estimate significance of progeny scores. Default value is 10000.
#'@param z_score if true, the z-scores will be returned for the pathway activity estimations. Else, the function returns a normalised z-score value between -1 and 1.
#'@param get_nulldist if true, the null score distribution used for normalisation will be returned along with the actual normalised score data frame.
#'@return This function returns a list of two elements. The first element is a dataframe of p*m+1 dimensions, where p is the number of progeny pathways, and m is the number of samples/contrasts.
#'Each cell represent the significance of a progeny pathway score for one sample/contrast. The signifcance ranges between -1 and 1. The significance is equal to x*2-1, x being the quantile of the progeny pathway score with respect to the null distribution.
#'Thus, this significance can be interpreted as the equivalent of 1-p.value (two sided test over an empirical distribution) with the sign indicating the direction of the regulation.
#'The sceond element is the null distribution list (a null distribution is generated for each sample/contrast).
runProgenyFast <- function(df,weight_matrix,k = 10000, z_scores = T, get_nulldist = F)
{
  resList <- list()
  if(get_nulldist)
  {
    nullDist_list <- list()
  }
  
    for(i in 2:length(df[1,]))
  {
    current_df <- df[,c(1,i)]
    current_df <- current_df[complete.cases(current_df),]
    t_values <- current_df[,2]
    
    current_weights <- weight_matrix
    
    names(current_df)[1] <- "ID"
    names(current_weights)[1] <- "ID"
    
    common_ids <- merge(current_df, current_weights, by = "ID")
    common_ids <- common_ids$ID
    common_ids <- as.character(common_ids)
    
    row.names(current_df) <- current_df$ID
    current_df <- as.data.frame(current_df[common_ids,-1])
    
    row.names(current_weights) <- current_weights$ID
    current_weights <- as.data.frame(current_weights[common_ids,-1])
    
    current_mat <- as.matrix(current_df)
    current_weights <- t(current_weights)
    
    scores <- as.data.frame(current_weights %*% current_mat)
    
    null_dist_t <- replicate(k, sample(t_values,length(current_mat[,1]), replace = F))
    
    null_dist_scores <- current_weights %*% null_dist_t
    
    if(get_nulldist)
    {
      nullDist_list[[i-1]] <- null_dist_scores
    }
    
    if(z_scores)
    {
      scores$mean <- apply(null_dist_scores,1,mean)
      scores$sd <- apply(null_dist_scores,1,sd)
      resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
    else
    {
      for(j in 1:length(weight_matrix[,-1]))
      {
        ecdf_function <- ecdf(null_dist_scores[j,])
        scores[j,1] <- ecdf_function(scores[j,1])
      }
      score_probas <- scores*2-1
      
      resListCurrent <- score_probas[,1]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
  }
  names(resList) <- names(df[,-1])
  resDf <- as.data.frame(resList)
  if(get_nulldist)
  {
    names(nullDist_list) <- names(df[,-1])
    return(list(resDf, nullDist_list))
  }
  else
  {
    return(resDf)
  }
}
#'@param df A data.frame of n*m+1 dimension, where n is the number of omic features to be considered and m is the number of samples/contrasts.
#'The first column should be the identifiers of the omic features. These identifiers must be coherent with the identifers of the weight matrix.
#'@param weight_matrix A progeny coeficient matrix. the first column should be the identifiers of the omic features, and should be coherent with the identifiers provided in df.
#'@param k The number of permutations to be preformed to generate the null-distribution used to estimate significance of progeny scores. Default value is 10000.
#'@param z_score if true, the z-scores will be returned for the pathway activity estimations. Else, the function returns a normalised z-score value between -1 and 1.
#'@param get_nulldist if true, the null score distribution used for normalisation will be returned along with the actual normalised score data frame.

#1. Prepare df:
progeny.df <- emat[common_genes, , drop = FALSE]
gene.names <- as.data.frame(row.names(progeny.df))
names(gene.names) <- "Gene"
progeny.df <- cbind(gene.names, progeny.df)
#2. Prepare weight matrix:
gene.names2 <- as.data.frame(row.names(model_matched))
names(gene.names2) <- "Gene"
progeny.cm <- cbind(gene.names2, model_matched)
#3. Apply permutations:
progeny.permutations <- t(runProgenyFast(progeny.df, progeny.cm, k = 10000, z_scores = T, get_nulldist = F))
##HEATMAP
pheatmap(progeny.permutations, scale = "column")
#Export heatmap as pdf
pdf(file = "figs/Heatmap_progeny.pdf", width = 5, height = 5)
pheatmap(progeny.permutations, scale = "column")
dev.off()

#Checking for differences between the groups:
#We check if the control is different to the treated condition using a linear model:
#1. Stablish control samples:
controls <- sample_info$cell_type == pop[1]
#2. Perform linear regression to calculate enrichment of pathways:
result1 <- apply(t(progeny.permutations), 1, function(x) {
  broom::tidy(lm(x ~ !controls)) %>%
    filter(term == "!controlsTRUE") %>%
    select(-term)
})
result1 <- mutate(bind_rows(result1), Pathway=names(result1))
results1 <- result1[,c(5,1:4)] 
write_xlsx(results1, "output/Progeny_permutations.xlsx")

#Plot pathways depending on significance and activation:
a <- as.data.frame(results1)
path.sig <- data.frame(stat = a$statistic, Sig = -log2(a$p.value))
row.names(path.sig) <- a$Pathway
Progeny_plot <- ggplot(path.sig, aes(x=stat, y=Sig)) +
  geom_point( aes(size=abs(stat),color=Sig),
              show.legend = F) +
  scale_color_gradient(low = "blue", high = "red")+
    geom_text_repel(data=subset(path.sig, Sig>6.6), 
                   label=rownames(subset(path.sig, Sig>6.6)),
                   size = 5,
                   show.legend = F) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 6.6) +
  labs(title=paste("PROGENy pathways activity", pop[1],"vs",pop[2]),
       x ="t value", y = "-log2(pvalue)") +
  theme(plot.title = element_text(face = "bold", colour = "black", size = 16),
        axis.title.x = element_text(color="navy", size=14, face="bold"),
        axis.title.y = element_text(color="navy", size=14, face="bold"))
#Export progeny dotplot as pdf
pdf(file = "figs/Dotplot_progeny.pdf", width = 5, height = 5)
Progeny_plot
dev.off()

#########################################################################################################
# @Aurelien's script: PROGENy can also be performed on a contrast (e.g. DESeq2 results, using log2FC)
#'\code{progenyScores}
#'#'This function compute progeny pathway scores as a series of weighted sum.
#'For each pathway a score is computed as the weighted sum of the sample/contrast statistic using the progeny coeficients as weights.
#'The gene identifiers between the measurments and the progeny coeficient matrix should be coherent.
#'#'@param df a n*m data frame, where n is the number of omic features (genes). m isn't really important, as long as at least one column corespond to a sample or contrast statistic. One of the column should correspond to the gene symboles.
#'@param cm a progeny coeficient matrix. One of the column should be the gene symboles.
#'@param dfIndex an integer corresponding to the column number of the gene identifiers of df.
#'@param FCIndex an integer corresponding to the column number that contains the statistic to be consdered to compute the progeny scores.
#'@param cmIndex an integer corresponding to the column number of the gene identifiers of the weight matrix.
#'@return a named vector where each element is a scores and the names are the corresponding pathways.
progenyScores <- function(df, cm, dfIndex = 1, FCIndex = 3, cmIndex = 1) {
  names(df)[dfIndex] <- "X1"
  names(cm)[cmIndex] <- "X1"
  df <- df[complete.cases(df[,FCIndex]),]
  merged <- merge(df[,c(dfIndex,FCIndex)],cm)
    for (pathway in names(cm[,-cmIndex]))
  {
    merged[,pathway] <- merged[,2]*merged[,pathway]
  }
  progeny_scores <- colSums(merged[,c(3:length(merged[1,]))])
  names(progeny_scores) <- names(merged[,c(3:length(merged[1,]))])
  return(progeny_scores)
}
# # Import Deseq2 results
# # if you want to generate the data here:
# # Remove genes with very low counts:
# dset1 <- dset[ rowSums (counts(dset)) > 10, ]
# # Set control population for the DE comparison:
# dset1$cell_type <- relevel(dset1$cell_type, pop[1])
# # Apply DESeq to do differential expression:
# dset1 <- DESeq(dset1)
# # Clean up genes with low counts and high variability with Shrinkage
# resLFC <- lfcShrink(dset1, coef=resultsNames(dset1)[2], type="apeglm")
# # Order genes based on their pvalue
# resOrdered <- as.data.frame(resLFC[order(resLFC$log2FoldChange, decreasing = T),])
# 
# write.table(resOrdered, paste0("data/DESeq2_",
#                           strsplit(resultsNames(dset1)[2], "_")[[1]][3], "_vs_",
#                           strsplit(resultsNames(dset1)[2], "_")[[1]][5], ".txt"),
#                           quote = F, dec = ".", sep = "\t")
# write_xlsx(resOrdered, paste0("data/DESeq2_",
#                           strsplit(resultsNames(dset1)[2], "_")[[1]][3], "_vs_",
#                           strsplit(resultsNames(dset1)[2], "_")[[1]][5], ".xlsx"))
 

gene.names <- as.data.frame(row.names(DESeq2))
progeny.df <- cbind(gene.names, DESeq2$log2FoldChange)
names(progeny.df) <- c("Gene", "log2FC")
progenyscores2 <- progenyScores(progeny.df, progeny.cm, dfIndex = 1, FCIndex = 2, cmIndex = 1)
PS2 <- as.data.frame(progenyscores2)
write_xlsx(PS2, "output/Progeny_scores_from_DESeq2.xlsx")
progeny2graph <- ggplot(PS2, aes(x=rownames(PS2), y=progenyscores2)) + geom_col(fill="navyblue") +
  xlab("Progeny pahtways") + ylab("Contrast score") +
  ggtitle(paste("Progeny analysis based on DEseq analysis",pop[1],"vs",pop[2])) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face="bold")) +
  theme(axis.text.x = element_text(angle = 60, size = 10, hjust =1, face="bold"))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.3, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"))

pdf(file = "figs/Graphbar_progenyon_DESeq.pdf", width = 8, height = 8)
progeny2graph
dev.off()
  
###########################################################################
#Scatter plots to check how genes within the pathways behave:
#LOOP FOR ALL PATHWAYS
Min0 = function(x) {
  min(x[x!=0])
}
pdf(file=paste("figs/Progeny score","pathways",".pdf"))
for (i in names(model)){
#Scatter plots to check how genes within the pathways behave:
#Retrieve genes that are involved in the pathway:
genes <- rownames(model[which(model[,i] != 0),])
#Filter the DESeq2 file to contain only these genes:
DESeq_ <- na.omit(DESeq2[genes,"log2FoldChange", drop = FALSE])
#Filter the progeny scoring matrix to contain only these proggenes too:
progeny_ <- na.omit(model[genes,i, drop = FALSE])
#Select common genes:
genes_common <- intersect(rownames(DESeq_), rownames(progeny_))
x = as.data.frame(DESeq_[genes_common,"log2FoldChange", drop = FALSE])
y = progeny_[genes_common,i, drop = FALSE]
stopifnot(rownames(x) == rownames(y))
.plot <- merge(x = x, y = y, by = "row.names")
row.names(.plot) <- .plot [,1]
.plot <- .plot[,2:3]
colnames(.plot) <- c("log2FC", "progeny")
#We need to use these ` ` symbols around the progeny weights variable
#because there is a space and this confuses R. Using these you mark it
#as a one single variable
# Change the point size depending on FC
#And + geom_point(aes(size=abs(log2FC)), show.legend = F)
#Add gene names to those that have a |FC|>2 and |progeny|>1, then change the size of the texts:
#And we add colour depending on FC, a gradient:
gplot <- ggplot(.plot, aes(x=log2FC, y=progeny)) + 
  geom_point(aes(size=abs(log2FC), color=abs(log2FC)), show.legend = F) +
  geom_label_repel(data=subset(.plot, abs(log2FC)>2 & abs(progeny)>1 ), 
                   label=rownames(subset(.plot, abs(log2FC) > 2 & abs(progeny)>1)),
                   size = 4, aes(fontface=2),
                   show.legend = FALSE)+
  scale_color_gradient(low = "blue", high = "red")+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(title=paste(i,"pathway",pop[1],"vs",pop[2]),
       x ="log2FoldChange", y = "PROGENy score",
       subtitle = paste0("Analysis p.value=",results[i,"p.value"], nsmall=3)) +
  theme(plot.title = element_text(face = "bold", colour = "black", size = 18),
        plot.subtitle = element_text(face = "bold", colour = "black", size = 14),
        axis.title.x = element_text(color="navyblue", size=14, face="bold"),
        axis.title.y = element_text(color="navyblue", size=14, face="bold"))+
        coord_cartesian(xlim = c(floor(Min0(DESeq2$log2FoldChange)), ceiling(max(DESeq2$log2FoldChange))),
                        ylim = c(floor(Min0(model)), ceiling(max(model))))
# If we wanted to change the fontface: 
#Allowed values : 1(normal), 2(bold), 3(italic), 4(bold.italic)
print(gplot)
}
dev.off()
####################################################################

