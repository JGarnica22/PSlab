#This script is to perform GSEA analysis with RNAseq data (can be used for other kind of data)
#It was created for R 3.6.2 version (2019-12-12)
#Copyright (C) 2020  Patricia Solé Sánchez and Mireia Ortega Crespo

# Check if required packages are installed, if not install:

bioc.packages <- c("fgsea", "qusage", "biomaRt", "stats")
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
library(stats)

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/GSEA")

## Specify parameters to be used along the script:
# Indicate name of txt file containing DESeq2 analysis result previously performed
Desq2file <- "DESeq2_TFH_vs_TH0.txt"
# Indicate populations of interest, to be compared:
# First indicate control population, then sample:
pop <- c("Th0", "TFH")
# Indicate species used in the analysis:
species <- "mouse"

# Indicate gene set file for the analysis:
gmt <- "h.all.v7.1.symbols.gmt"


##You need to work with differential expression data from your previous analysis.
# Read the DESeq2 results:
DESeq2 <- read.table (file = paste0("data/",Desq2file),
                      sep = "\t", 
                      quote = "",
                      dec = ".", 
                      header = T, row.names = "Gene")

# Read gene sets (gmt.file):
#The gmt.file is a list object that contains a set of vectors corresponding to the
#gene sets. Each geneset (vector), will contain the genes within the geneset.
gmt.file <- read.gmt(paste0("data/", gmt))
gmt.file

# First, we need to create a ranking of genes for GSEA. 
#You can rank genes based on FC or on significance (log10 of the padj), we will use preferably FC

#1) Order genes depending on FC:
resOrdered <- DESeq2[order(DESeq2$log2FoldChange, decreasing = TRUE),]
resOrdered

#2) Ranks object is a vector of FC with the names of the genes:
ranks <- resOrdered[, "log2FoldChange"]
names(ranks) <- row.names(resOrdered)


#Gene names in the gmt file are human genes (capital letters)
#If you are using mouse data (specified at the beginning of the script),
#this next code will transform the gmt list to mouse genes:

#Create the BioMart:
if (species == "mouse"){
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
  #Use for loop to change the names of all the genes in all the gene sets (for within for):
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
  gmt.file <- gmt.mouse
}

# Run fGSEA:
fgseaRes <- fgsea(pathways = gmt.file, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
fgseaRes

#Table plots for all hallmarks:
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf(file = "figs/Table_GSEA.pdf", width = 11, height = 5)
plotGseaTable(gmt.mouse[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#If from the plot above one can see that there are very similar pathways in the table (for example
#Mitotic_Metaphase_and_Anaphase and Mitotic_Anaphase). To select only independent pathways one 
#can use collapsePathways function:
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      gmt.mouse, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

pdf(file = "figs/Table_GSEA_collapsed.pdf", width = 11, height = 5)
plotGseaTable(gmt.mouse[mainPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#To save the results data:table::fwrite function can be used:
fwrite(fgseaRes, 
       file= paste0("output/fgseaRes_", pop[2], "_vs_", pop[1], ".txt"),
       sep="\t", sep2=c("", " ", ""))

#Plot GSEA for an specific gene set:
#Indicate gene set of interest
geneset <- "HALLMARK_G2M_CHECKPOINT"

pdf(file = paste0("figs/GSEA_", geneset, ".pdf"), width = 5, height = 5)
plotEnrichment(gmt.mouse[[geneset]],
               ranks) + 
  labs(title="G2M_CHECKPOINT")
dev.off()


#Plot GSEA for all gene sets in a single pdf:
pdf(paste0("figs/GSEA_all_plots.pdf"), height=5, width=6)
for (i in 1:length(gmt.file)){
  plt <- plotEnrichment(gmt.mouse[[i]],
                        stats= ranks,
                        gseaParam = 1, 
                        ticksSize = 0.5) + 
    labs(title= names(gmt.mouse[i])) + 
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  print(plt)
}
dev.off()
