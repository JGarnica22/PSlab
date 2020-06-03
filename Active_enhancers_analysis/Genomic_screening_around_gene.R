#This script is to look for annotated regions around genes of interest
#It was created for R 3.6.3 version (2020-05-22)
#Copyright (C) 2020  Patricia Solé Sánchez and Josep Garnica Caparrós
#####################################################################
# Check if required packages are installed, if not install:
cran.packages <- c("tidyr", "stringr", "Cairo", "snakecase", 
                   "plyr", "dplyr", "writexl",
                   "ggplot2", "gridExtra")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}
bioc.packages <- c("GenomicRanges", "trackViewer", "biomaRt",
                   "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db",
                   "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

library(GenomicRanges)
library(trackViewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyr)
library(stringr)
library(Cairo)
library(snakecase)
library(plyr)
library(dplyr)
library(biomaRt)
library(writexl)
library(ggplot2)
library(gridExtra)

# Set your working directory (the project you are working in):
setwd("/home/jgarnica/R/GenomicRanges_Active_enhancers")
setwd("/Users/patri/Documents/LAB/TESIS DOCTORAL 2015-2020/INTERNSHIP/7_GenomicRanges")

# Indicate populations to work with
pop <- c("Tconv", "TR1")

# Indicate species working with (mouse or human)
species <- "mouse"

#Generate database for species to be studied
if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
}

#This analysis looks for certain elements (such as active enhancers) around 100kb from genes of interest.
#Set the genes to study or import file:
genes <- c("Il10", "Il21")

# Create a empty list:
Act.enh.gene <- vector(mode = "list",length = length(genes))
# Set the names of the gene sets in the list:
names(Act.enh.gene) <- genes

for (g in names(Act.enh.gene)) {
  enhancers <- vector()
  id <- get(g, org.SYMBOL2EG)
  gr <- genes(TxDb)[id]
  #set the window width of your search, as a defect is set at 100.000 bp
  gr100kb <- resize(gr, width(gr)+100000, fix = "center")
  
  #Introduce here your gr file (or generated such as the example for active enhancers
  #overlaping ATAC and ChIP) annotating the position of the elements you want to search around your genes
  #act.enhan <- file_of_reference
  #subsetByOverlaps(eval(as.symbol(grep(paste0("ChIP.",pop[i]), names(.GlobalEnv),value=TRUE))), 
  #eval(as.symbol(grep(paste0("ATAC.",pop[i]), names(.GlobalEnv),value=TRUE))))
  #Filter overlapping peaks in promoters:
  #prom <- promoters(TxDb)
  #act.enhan <- act.enhan[-c(unique(queryHits(findOverlaps(act.enhan, prom))))]
  act.enh.by.gene <- subsetByOverlaps(act.enhan, gr100kb)
  if (length(act.enh.by.gene)>0){
    df <- data.frame(gene = rep(g, length(act.enh.by.gene)),
                     location = paste0(seqnames(act.enh.by.gene),":", ranges(act.enh.by.gene)))
    df <- df[order(df$location),]
    df <- ddply(df, .(gene), summarize, 
                gene=paste(unique(gene),collapse=","),
                location=paste(unique(location),collapse=","))
    enhancers <- df[,"location"]
    Act.enh.gene[[g]] <- enhancers
  }
  
  if (file.exists("Act.enh.gene")){
    trial <-as.data.frame(unlist(Act.enh.gene))
    write.table(trial, paste0("output/", "your_element" ,"_around_genes.txt"),
                sep = "\t", row.names = T, col.names = F, quote = F)
    #rm(Act.enhancers)
  }
}
