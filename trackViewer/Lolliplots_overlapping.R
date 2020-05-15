#This script is to visualize methylome data by trackViewer doing lolliplots
#It was created for R 3.6.3 version (2020-05-03)
#Copyright (C) 2020  Patricia Sole Sanchez
##########################################################
# Check if required packages are installed, if not install:
cran.packages <- c("Cairo", "stringr", "tidyr")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}
bioc.packages <- c("GenomicRanges", "Gviz", "trackViewer", "rtracklayer", 
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
library(Gviz)
library(trackViewer)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyr)
library(stringr)
library(Cairo)
library(plotrix)

# Set your working directory (the project you are working in):
setwd("/home/jgarnica/R/trackviewer")

# Indicate populations to work with, first indicate control population.
pop <- c("Tconv", "Tet")

#Inidicate control and sample names
comp <- c("Ctl", "Smp")

# Indicate species working with (mouse or human)
species <- "mouse"

# Indicate the genes to visualize
plots <- c("Il10", "Il21", "Irf4", "Pdcd1")
# or read a txt file containing a list with all the plots

if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
}

#Import file with methylation information for each nucleotid comparing the two populations to compare
file_list_txt <- list.files(path="/home/jgarnica/R/trackviewer/data",pattern= "CG.")

for (file in c(file_list_txt)){
  x <- read.table(paste0("data/", file), sep = "\t", dec = ".", header = T, quote = "")
  for (u in 1:length(comp)){
  y <- GRanges(seqnames = x$chr, 
               ranges = IRanges(start = x$pos,
                                width = 1), 
               strand = NULL,
               #We cross here columns because we assume groups were crossed
               #in the analysis!!
               score = as.integer(x[,u+2]*100))
  assign(paste0("lolli.", comp[u]), y)
    }
  

#Lolliplots overlapping samples:
#Call all GRanges objects to plot
grl <- c(mget(grep("lolli*",names(.GlobalEnv),value=TRUE)))
#set colors to use
colors <- rainbow(length(pop))
pdf(paste0("figs/Lolliplot_", paste0(pop, collapse ="_"), ".pdf"), width = 12, height = 6)
for (plot in plots) {
  id <- get(plot, org.SYMBOL2EG)
  gr <- genes(TxDb)[id]
  gene <- geneTrack(id,TxDb)[[1]]$dat
  gene$fill <- "lightblue"
  gene$height <- 0.03
  grlo <- GRanges()
  for (p in 1:length(grl)){
    grl[[p]]$color <- colors[p]
    grl[[p]]$border <- "gray30"
    grlo <- c(grlo, grl[[p]])
  }
  xaxis <- c(seq(from = start(gr), 
                 to = end(gr), 
                 by = 500))
  yaxis <- c(0, 100)
  legends <- list(list(labels=pop, 
                       fill= colors,
                       col = c("gray30","gray30")))
  
  lolliplot(grlo, gene, type = "circle", 
            xaxis= xaxis,
            yaxis= yaxis,
            legend= legends,
            ylab= plot)
    }
dev.off()
}
