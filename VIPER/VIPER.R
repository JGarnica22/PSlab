#This script is to infer TF activity in our test samples running msVIPER function
#It was created for R 3.6.0 version (2019-04-26)
#Copyright (C) 2020  Patricia Solé Sánchez and Mireia Ortega Crespo

## APPLY VIPER TO INFER PROTEIN ACTIVITY BASED ON OUR ARACNe-DERIVED REGULATORY NETWORK 
# VIPER compute TF activity changes from a differential gene expression signature:

#Check if required packages are installed, if not install:

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

# Indicate name of txt file containing rlog counts from the DESeq2 analysis:
rlog.norm.counts <- "rlog.norm.counts.txt"

# Indicate name of the txt file containing normalized matrix for CD4.
# Remember that tissue-specific data needs to be downloaded from GEO
norm.matrix.CD4 <- "norm.matrix.CD4.txt"

# Indicate name of txt file containing ARACNe output .adj file (the expression dataset used 
#by ARACNe to reverse engineer the network):
adjfile <- "network.txt"

# Indicate populations of interest, to be compared:
# First indicate control population, then sample:
pop <- c("Th0", "TFH")

# Read differential expression data from DESeq2:
DESeq2 <- read.table (file = paste0("data/",Desq2file),
                           sep = "\t", 
                           quote = "",
                           dec = ".", 
                           header = T, row.names = NULL)

# Exclude probes with unknown or duplicated gene symbol:
DESeq2 = subset(DESeq2, Gene != "" )
DESeq2 = subset(DESeq2, ! duplicated(Gene))

# Read expression data used for ARACNe:
exprdata <- read.table(file = paste0("data/",norm.matrix.CD4),
                       header = T, sep = "\t", dec = ".", quote = "",
                       row.names = NULL, stringsAsFactors = FALSE)

# Read rlog.norm.counts file from DESeq2:
rlog.norm.counts2 <- read.table(paste0("data/",rlog.norm.counts), 
                                header = T, sep = "\t", dec = ".", quote = "")

# Convert "Gene" column to row.names:
rlog.norm.counts3 <- data.frame(rlog.norm.counts2, row.names = "Gene")


# Estimate z-score values for the GES (Gene Expression Signature):
myStatistics = matrix(DESeq2$log2FoldChange, dimnames = list(DESeq2$Gene, 'logFC') )
myPvalue = matrix(DESeq2$pvalue, dimnames = list(DESeq2$Gene, 'pvalue') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]

# Generate the regulon object from ARACNe network:
rownames(exprdata) <- exprdata$gene
exprdata <- exprdata[,colnames(exprdata)!="gene"]
exprdata <- as.matrix(exprdata)
class(exprdata)

# Convert dset to class Expression set:
dset <- ExpressionSet(assayData=exprdata)
class(dset)

# Generate the regulon to be used in VIPER:
regulon <- aracne2regulon(paste0("data/",adjfile), dset, gene = FALSE, verbose = TRUE)
print(regulon)

# Create an histogram to check regulon distribution:
pdf(file = "figs/Histogram.pdf", width = 11, height = 5)
hist(unlist(lapply(regulon, function(TF) length(TF$tfmode))),
     main="Histogram", xlab="No.targets", breaks=25,las=1)
dev.off()

# Create a matrix from rlog.norm.counts to perform null distribution:
MAT <- as.matrix(rlog.norm.counts3)
DF <- data.frame(sample = colnames(MAT), group = "", stringsAsFactors = FALSE)
for(i_group in pop){
  DF[grepl(i_group,DF$sample), "group"] <- i_group
}

row.names(DF) <- DF[,1]
DF <- DF[,2,  drop = FALSE]

eset <-  ExpressionSet(assayData=MAT, phenoData=new("AnnotatedDataFrame",data=DF))

#Run nullmodel with ttestNull function:
nullmodel <- ttestNull(eset, "group", pop[2], pop[1], per = 1000,
                       repos = TRUE, verbose = TRUE)

# Estimate TF activities with msviper:
mrs <- msviper(mySignature, regulon, nullmodel, verbose = TRUE)
summary(mrs)

TFact = data.frame(Regulon = names(mrs$es$nes),
                   Size = mrs$es$size[ names(mrs$es$nes) ], 
                   NES = mrs$es$nes, 
                   p.value = mrs$es$p.value, 
                   FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TFact = TFact[ order(TFact$p.value), ]

#Export results of TF activity as a txt file:
write.table(TFact, "output/TFactivity_test.txt",
            sep = "\t", dec = ".",
            quote = F, row.names = F, col.names = T)

# Plot a Volcano plot FALTARIA POSAR MÉS MACO:
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
