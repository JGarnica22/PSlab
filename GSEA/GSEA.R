#This script is to perform GSEA analysis with RNAseq data (can be used fro other kind of data)
#It was created for R 3.6.2 version (2019-12-12)
#Copyright (C) 2020  Patricia Solé Sánchez and Mireia Ortega Crespo

# Check if required packages are installed, if not install:

bioc.packages <- c("fgsea", "qusage", "biomaRt")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

cran.packages <- c("ggplot2", "data.table")
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
library(fgsea)
library(qusage)
library(biomaRt)
library(ggplot2)
library(data.table)

# Set your working directory (the project you are working in):

setwd("/Users/Mireia/Desktop/R_Home")

#You need to work with differential expression data from your previous analysis.
#First, we need to create a ranking of genes for GSEA. 
#Then, we will order genes based on FC and also on significance (log10 of the padj):

#1) Order genes depending on FC:
DGE.LFC <- read.table("/Users/Mireia/Desktop/R_Home/Differential_expression_examen2.txt")
DGE.FC.ord <- DGE.LFC[order(DGE.LFC$log2FoldChange, decreasing = TRUE),]
DGE.FC.ord

#2) Ranks object is a vector of FC with the names of the genes
ranks <- DGE.FC.ord[, "log2FoldChange"]
names(ranks) <- row.names(DGE.FC.ord)

#Read hallmarks downloaded from MSigDB
#The gmt.file is a list object that contains a set of vectors corresponding to the
#gene sets. Each geneset (vector), will contain the genes within the geneset.

gmt.file <- read.gmt("hallmarks.gmt")
gmt.file

#Read a Hallmak as an exaple: HALLMARK_G2M_CHECKPOINT 
gmt.file["HALLMARK_G2M_CHECKPOINT"]

#Gene names in the gmt file are in capital letters (human genes), while gene names in "ranks" 
#are not (mouse genes). 
#The CORRECT way is using Biomart to find orthologs for the gene names in the 
#gene sets:
#Create the BioMart
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

BM1 <- getBM(attributes=c("ensembl_gene_id","gene_biotype",
                          "external_gene_name"), mart = mouse, verbose = T)
BM2 <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), 
             mart = mouse, verbose = T)
#Merge two data frames in one by the common column ("ensembl_gene_id")
BM <- merge(x = BM1, y = BM2, by = "ensembl_gene_id", all.x = T)
BM <- unique(BM)
BM <- BM[,c(3:4)]
names(BM) <- c("mouse", "human")

#Use the "for" function or loop to change the names of all the genes in all the gene sets (for within for):

# Create a empty list:
gmt.mouse <- vector(mode = "list",length = length(gmt.file))

# Set the names of the gene sets in the list:
names(gmt.mouse) <- names(gmt.file) 

#Change the gene names from human to mouse: 
#For each Gene_set we create an empty vector, that will after contain the genes
for(Gene_set in names(gmt.file)) {
  geneset_mouse <- vector()
  #Then for each gene in each gene set, we create a vector with the mouse gene names and we move this data to the
  #vector we created on top for each gene set (named geneset_mouse)
  for (h_gene in gmt.file[[Gene_set]]){
    m_gene <- BM[which(BM$human==h_gene),"mouse"]
    geneset_mouse <- c(geneset_mouse, m_gene)
  }
  #And finally we need to copy the information in each vector to the empty list - with names- list we created before (gmt.mouse)
  #We replace each element of the list (that was empty, NULL) with the vector we created with the mouse gene names
  gmt.mouse[[Gene_set]] <- geneset_mouse
}

#Run fGSEA:
fgseaRes <- fgsea(pathways = gmt.mouse, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
fgseaRes

#Table plots for all hallmarks
pdf(file = "figs/table_GSEA.pdf", width = 11, height = 5)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gmt.mouse[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#If from the plot above one can see that there are very similar pathways in the table (for example
#Mitotic_Metaphase_and_Anaphase  and Mitotic_Anaphase). To select only independent pathways one 
#can use collapsePathways function:

pdf(file = "figs/table_GSEA_collapsed.pdf", width = 11, height = 5)
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      gmt.mouse, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(gmt.mouse[mainPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#To save the results data:table::fwrite function can be used:

fwrite(fgseaRes, file="output/fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

#Plot GSEA for an specific hallmark
pdf(file = "figs/GSEA.pdf", width = 5, height = 5)
plotEnrichment(gmt.mouse[["HALLMARK_G2M_CHECKPOINT"]],
               ranks) + 
  labs(title="G2M_CHECKPOINT")
dev.off()


###FALTA###
#Plot GSEA for all hallmarks in a single pdf
for (i in gmt.mouse){
  pdf(paste0(i,"figs/GSEA_all.pdf"),height=5,width=5) 
  plt <- plotEnrichment(pathway = ranks[[i]], 
                        gseaParam = 1, ticksSize = 0.5, stats= ranks) + 
    labs(title=i) + theme(plot.title = element_text(hjust = 0.5, face="bold"))
  print(plt)
  dev.off()
}

