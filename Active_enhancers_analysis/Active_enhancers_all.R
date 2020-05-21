#This script is to infer active enhancer regions based on ATAC anc ChIP-seq data
#It was created for R 3.6.3 version (2020-05-22)
#Copyright (C) 2020  Patricia Solé Sánchez and Josep Garnica Caparrós
#####################################################################
# Check if required packages are installed, if not install:
cran.packages <- c("Cairo", "stringr", "tidyr", "snakecase", "plyr")
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

# Set your working directory (the project you are working in):
setwd("/home/jgarnica/R/GenomicRanges_Active_enhancers")

## NOTE: CHANGE NAMES OF THE FILES TO AVOID CONFUSIONS, IN THE DATA FILES SHOULD APPEAR THE NAME OF THE TECHNIQUE AND ONLY THE POPULATION STUDIED

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

#Load all files needed
#Load DMR file between two samples:
DMR <- read.table(paste0("data/", list.files(path="/home/jgarnica/R/GenomicRanges_Active_enhancers/data", pattern= "DMR")),
                  sep = "\t", quote = "",
                  dec = ".", header = T, na.strings = T)

#Load and prepare sharted OCR between two populations:
socr <- read.table(paste0("data/", list.files(path="/home/jgarnica/R/GenomicRanges_Active_enhancers/data", pattern= "OCR")),
                   sep = "\t", dec = ".",header = TRUE, quote = "", stringsAsFactors = F)
socr <- socr[, "Region.ID", drop = F]
socr$Chr <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[1]))
socr$ranges <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[2]))
grocr <- GRanges(seqnames = socr$Chr, 
               ranges = socr$ranges, 
               strand = NULL)

for (i in c(1:length(pop))) {
  for (m in c("ChIP", "ATAC")){
  file_list <- list.files(path="/home/jgarnica/R/GenomicRanges_Active_enhancers/data", pattern= m)
  #Read file tables
  tble <- read.table(paste0("data/", grep(pop[i], file_list, value = T)),sep = "\t", quote = "",
                   dec = ".", header = T, na.strings = T)
  if (m == "ChIP"){
  tble[, c("Sample.name", "Absolute.summit", "Pileup", "X.log10.qvalue.",
           "Peak.name", "Transcript.IDs")] <- NULL
  tble[, 14:ncol(tble)] <- NULL
  
  names(tble) <- c("Chr", "Start", "End", "Length", "-log10.pval", "FoldEnrichment",
                   "Anno.Gene", "Strand", "Transcript.start", "Transcript.end",
                   "Gene.section", "Distance.to.TSS", "Distance.to.TTS")
  tble <- tble[order(as.numeric(gsub("chr", "", tble$Chr)), 
                     as.numeric(tble$Start),
                     decreasing = F, na.last = T), ]
  
  gr <- GRanges(seqnames = tble$Chr, 
                ranges = paste0(tble$Start,"-",tble$End), 
                strand = NULL,
                `-log10.pval`= tble$`-log10.pval`,
                FoldEnrichment = tble$FoldEnrichment,
                Anno.Gene = tble$Anno.Gene,
                Peak.location = tble$Gene.section,
                Distance.to.TSS = tble$Distance.to.TSS)
  grchip <- GRanges(seqnames = tble$Chr, 
                ranges = paste0(tble$Start,"-",tble$End), 
                strand = NULL,
                `-log10.pval`= tble$`-log10.pval`,
                FoldEnrichment = tble$FoldEnrichment,
                Anno.Gene = tble$Anno.Gene,
                Peak.location = tble$Gene.section,
                Distance.to.TSS = tble$Distance.to.TSS)
  } else {
    tble <- tble[, c(2:4,12)]
    names(tble) <- c("Chr", "Start", "End", "Anno.Gene")
    atac <- tble[order(as.numeric(gsub("chr", "", tble$Chr)), 
                       as.numeric(tble$Start),
                       decreasing = F, na.last = T), ]
    
    gr <- GRanges(seqnames = atac$Chr, 
                   ranges = paste0(atac$Start,"-",atac$End), 
                   strand = NULL,
                   Anno.Gene = atac$Anno.Gene)
  }

  assign(paste0(m,".", pop[i], ".gr"), gr)
  }
  
  
  #Find overlapping peaks
  #Careful:apparently order of objects in `findOverlaps` matters!
  overlap <- findOverlaps(eval(as.symbol(grep(paste0("ChIP.",pop[i]), names(.GlobalEnv),value=TRUE))), 
                          eval(as.symbol(grep(paste0("ATAC.",pop[i]), names(.GlobalEnv),value=TRUE))))
  olpeaks <- atac[unique(subjectHits(overlap)),]
  gr3 <- GRanges(seqnames = olpeaks$Chr, 
                 ranges = paste0(olpeaks$Start,"-",olpeaks$End), 
                 strand = NULL)
  write.table(olpeaks, file = paste0("output/", pop[i], "_ATAC_Overlapping_peaks_with_H3K27ac_ChIP.txt"),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  #Obtain active enhancers by filtering overlapping peaks in promoters:
  prom <- promoters(TxDb)
  inpromoters <- findOverlaps(prom, gr3)
  act.enh <- olpeaks[-c(unique(subjectHits(inpromoters))), 1:3]
  
  #Export files in desired formats
  for (o in c(".txt", ".bed")){
  write.table(act.enh, file = paste0("output/", pop[i] ,"_Active_enhancers", o),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  #Methylation in active enhancers
  gr4 <- GRanges(seqnames = act.enh$Chr, 
                 ranges = paste0(act.enh$Start,"-",act.enh$End), 
                 strand = NULL)
  names(DMR) <- c("Chr", "Start", "End", pop[1], pop[2])
  DMR <- DMR[order(as.numeric(gsub("chr", "", DMR$Chr)), 
                   as.numeric(DMR$Start),
                   decreasing = F, na.last = T), ]
  
  gr5 <- GRanges(seqnames = DMR$Chr, 
                 ranges = paste0(DMR$Start,"-",DMR$End), 
                 strand = NULL,
                 Met.smp = DMR[,pop[2]],
                 Met.ctl = DMR[,pop[1]])
  overlap <- findOverlaps(gr5, gr4)
  act.enh.DMR <- act.enh[unique(subjectHits(overlap)),]
  write.table(act.enh.DMR, file = paste0("output/", pop[i] ,"_Active_enhancers_with_DMR", o),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  #do the overlap in the other direction
  overlap2 <- findOverlaps(gr4, gr5)
  DMR.act.enh <- DMR[unique(subjectHits(overlap2)),]
  write.table(DMR.act.enh, paste0("output/", pop[i] ,"_DMR_Overlapping_Active_enhancers", o),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  
  #Methylation in H3k27ac
  #Find overlapping peaks (population-specific H3K27ac mark + open region)
  overlap3 <- findOverlaps(grchip, grocr)
  openH3K27ac <- socr[unique(subjectHits(overlap3)),]
  grH3 <- GRanges(seqnames = openH3K27ac$Chr, 
                 ranges = openH3K27ac$ranges, 
                 strand = NULL)
  openH3K27ac <- separate(openH3K27ac, col = "ranges", into = c("Start", "End"), sep = "-", remove = T)
  openH3K27ac <- openH3K27ac[, c("Chr", "Start", "End")]
  write.table(openH3K27ac, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac", o),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  
  #Obtain active enhancers by filtering overlapping peaks in promoters:
  inpromoters <- findOverlaps(prom, grH3)
  openH3K27ac <- openH3K27ac[-c(unique(subjectHits(inpromoters))), 1:3]
  write.table(openH3K27ac, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac_not_promoter", o),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  overlap5 <- findOverlaps(gr5, grH3)
  H3K27open.DMR <- openH3K27ac[unique(subjectHits(overlap5)),]
  H3K27open.DMR <- H3K27open.DMR[which(H3K27open.DMR$Chr!="NA"),]
  write.table(H3K27open.DMR, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac_not_promoter_with_DMR", o),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  
  overlap6 <- findOverlaps(grH3, gr5)
  DMR.H3K27open <- DMR[unique(subjectHits(overlap6)),]
  write.table(DMR.H3K27open, file = paste0("output/", pop[i] ,"_DMR_Overlapping_shared_ATAC_H3K27ac_not_promoter", o),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  }
}

