#This script is to infer TF activity in our test samples running msVIPER
#It was created for R 3.6.0 version (2019-04-26)
#Copyright (C) 2020  Patricia Solé Sánchez and Mireia Ortega Crespo

## APPLY VIPER TO INFER PROTEIN ACTIVITY BASED ON OUR ARACNe-DERIVED REGULATORY NETWORK 
# VIPER compute TF activity changes from a differential gene expression signature:

# Check if required packages are installed, if not install:

bioc.packages <- c("viper", "ggplot2")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

cran.packages <- c("dplyr")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}

# Load packages:
library(viper)
library(ggplot2)
library(dplyr)

# Set your working directory (the project you are working in):
setwd("/Users/Mireia/Desktop/R_Home/VIPER2")

## Specify parameters to be used along the script:

# Indicate name of txt file containing DESeq2 analysis result previously performed:
Desq2file <- "DESeq2_TFH_vs_TH0.txt"

# Indicate name of txt file containing rlog counts from teh DESeq2 script:
rlog_counts_test <- "rlog_normalized_DESeq.txt"
rlog_counts_test2 <- read.table(paste0("data/",rlog_counts_test), header = T, sep = "\t", dec = ".", quote = "")
#Convert Gene column to row.names:
rlog_counts_test3 <- data.frame(rlog_counts_test2, row.names = "Gene")
#Eliminate last conditions but this should not appear?¿?¿
rlog_counts_test4 <- dplyr::select(rlog_counts_test3, -c("Tet1", "Tet2", "Tet3", "Tet4"))

# Indicate name of txt file containing normalized matrix for CD4 (tissue-specific data needs to be downloaded from GEO):
#We can use ARSCH4 web to download expression data from mouse CD4+ T cell samples.
norm.matrix.CD4 <- "norm.matrix.CD4.txt"

# Indicate name of txt file containing ARACNe output .adj file (and the expression dataset used by ARACNe to reverse engineer the network):
adjfile <- "network.txt"

# Indicate populations of interest, to be compared:
# First indicate control population, then sample:
pop <- c("Th0", "TFH")

# Read differential expression data from DESeq2:
DEsignature <- read.table (file = paste0("data/",Desq2file),
                           sep = "\t", 
                           quote = "",
                           dec = ".", 
                           header = T, row.names = NULL)

# Exclude probes with unknown or duplicated gene symbol:
DEsignature = subset(DEsignature, Gene != "" )
DEsignature = subset(DEsignature, ! duplicated(Gene))

# Estimate z-score values for the GES (Gene Expression Signature):
#Check VIPER manual for details.

myStatistics = matrix(DEsignature$log2FoldChange, dimnames = list(DEsignature$Gene, 'logFC') )
myPvalue = matrix(DEsignature$pvalue, dimnames = list(DEsignature$Gene, 'pvalue') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]

#Generate the regulon object from the ARACNe network:
#Regulon objects can be generated from networks reverse engineered with the ARACNe algorithm. 
#This is performed by the function ´aracne2regulon´, which takes two arguments as input: the
#ARACNe output .adj file, and the expression dataset used by ARACNe to reverse engineer the network

# Read expression data used for ARACNe:
exprdata <- read.table(file = paste0("data/",norm.matrix.CD4),
                       header = T, sep = "\t", dec = ".", quote = "",
                       row.names = NULL, stringsAsFactors = FALSE)
rownames(exprdata) <- exprdata$gene
exprdata <- exprdata[,colnames(exprdata)!="gene"]
exprdata <- as.matrix(exprdata)
class(exprdata)

# Convert dset to class Expression set:
dset <- ExpressionSet(assayData=exprdata)
class(dset)

# Generate regulon to use in VIPER:
regulon <- aracne2regulon(paste0("data/",adjfile), dset, gene = FALSE, verbose = TRUE)
print(regulon)

# Check regulon distribution and export histogram:
pdf(file = "figs/Histogram.pdf", width = 11, height = 5)
hist(unlist(lapply(regulon, function(TF) length(TF$tfmode))),
     main="Histogram", xlab="No.targets", breaks=25,las=1)
dev.off()

# Create a null distribution:
DF <- data.frame(sample = colnames(MAT), group = "", stringsAsFactors = FALSE)
for(i_group in pop){
  DF[grepl(i_group,DF$sample), "group"] <- i_group
}

# Create phenoData:
row.names(DF) <- DF[,1]
DF <- DF[,2,  drop = FALSE]
eset <-  ExpressionSet(assayData=MAT, phenoData=new("AnnotatedDataFrame",data=DF))

#Where MAT is a matrix with gene expression values, rows=genes and cols=samples and DF
#is a data.frame with rows = samples, and cols=metadata for samples. For instance, group.
#Both cols from MAT and rows from DF must be the same and follow same order

nullmodel <- ttestNull(eset, "group", pop[2], pop[1], per = 1000,
                       repos = TRUE, verbose = TRUE)


# Estimate TF activities:
mrs <- msviper(mySignature, regulon, nullmodel, verbose = TRUE)
summary(mrs)

TFact = data.frame(Regulon = names(mrs$es$nes),
                   Size = mrs$es$size[ names(mrs$es$nes) ], 
                   NES = mrs$es$nes, 
                   p.value = mrs$es$p.value, 
                   FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TFact = TFact[ order(TFact$p.value), ]

write.table(TFact, "output/TFactivity_test.txt",
            sep = "\t", dec = ".",
            quote = F, row.names = F, col.names = T)

# Plot a Volcano plot FALTA FROMAT:
pdf(file = "figs/Volcano.pdf", width = 11, height = 5)
vol <- TFact[, c("NES", "FDR")]
vol[,2] <- -log10(vol[,2])
names(vol) <- c("NES", "-log10(FDR)")
Plot <- ggplot(vol, aes(x=NES, y=`-log10(FDR)`))
Plot + geom_point( aes(size=abs(NES),
                       color=`-log10(FDR)` > 1.3),
                   show.legend = F) + 
  scale_color_manual(values = c("TRUE" = "red","FALSE"="black")) +
  geom_vline(xintercept = c(-2, 2)) +
  geom_hline(yintercept = 1.3) +
  labs(title="TF inferred activities \n based on ARACNe network") +
  theme(plot.title = element_text(face = "bold", colour = "black", size = 14),
        axis.title.x = element_text(color="#993333", size=12, face="bold"),
        axis.title.y = element_text(color="#993333", size=12, face="bold"))
dev.off()
