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

# Indicate species working with (mouse or human)
species <- "mouse"

# Generate database for species to be studied
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

# Set the width of your search in base pairs
wi <- 100000 


# Call files containing the elements you want to look for (around genes) in .txt format
file.ele <- "TR1_Active_enhancers1"

tble <- read.table(paste0("output/", file.ele, ".txt"),sep = "\t", quote = "",
                   dec = ".", header = T, na.strings = T)
# Convert table into GRanges object
chr.loc <- GRanges(seqnames = tble$Chr, 
                   ranges = paste0(tble$Start,"-",tble$End), 
                   strand = NULL)

# Create a empty dataframes:
df_mr <- data.frame(matrix(ncol = 2, nrow = 0))
names(df_mr) <- c("gene","location")
df_su <- data.frame(matrix(ncol = 2, nrow = 0))
names(df_su) <- c("gene","location")

for (g in genes) {
  id <- get(g, org.SYMBOL2EG)
  gr <- genes(TxDb)[id]
  #set the window width of your search, as a defect is set at 100.000 bp
  gr100kb <- resize(gr, width(gr)+wi, fix = "center")
  chr.loc.by.gene <- subsetByOverlaps(chr.loc, gr100kb)
  if (length(chr.loc.by.gene)>0){
    df <- data.frame(gene = rep(g, length(chr.loc.by.gene)),
                     location = paste0(seqnames(chr.loc.by.gene),":", ranges(chr.loc.by.gene)))
    df <- df[order(df$location),]
    df2 <- ddply(df, .(gene), summarize, 
                 gene=paste(unique(gene),collapse=","),
                 location=paste(unique(location),collapse=","))
  } else {
    df <- data.frame(gene = g, location = paste0("No elements found ", wi, " bp around"))
    df2 <- data.frame(gene = g, location = paste0("No elements found ", wi, " bp around"))
  }
  # data frame including one row per each element found
  df_mr <- rbind(df_mr, df)
  # data frame with one row per gene and all elements listed in the same cell
  df_su <- rbind(df_su, df2)
  
}

write.table(df_mr, paste0("output/", file.ele ,"_around_genes_more_rows.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)
write.table(df_su, paste0("output/", file.ele ,"_around_genes_collapsed.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)
