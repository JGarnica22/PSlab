#This script is to visualize genomic data by trackviewer
#It was created for R 3.6.3 version (2020-05-03)
#Copyright (C) 2020  Patricia Sole Sanchez
##########################################################
# Check if required packages are installed, if not install:
cran.packages <- c("ggplot2", "ggrepel", "writexl", "pheatmap", "Cairo", "stringr", "tidyr")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}
bioc.packages <- c("DESeq2", "biomaRt", "trackViewer", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db", "GenomicRanges", "Gviz")
for (i in bioc.packages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i)
    print(paste(i,"just installed"))
  } else {
    print(paste(i,"was already installed"))
  }
}

library(GenomicRanges)
library(biomaRt)
library(trackViewer)
library(Gviz)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(tidyr)
library(stringr)
library(Cairo)
#add packages and respective libraries for human!!

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/trackviewer")

#Indicate populations to work with
pop <- c("Tconv", "Tet")

#Indicate species working with (typicall mouse or human)
species <- "mouse"
# Indicate the genes to visualize, for mouse minus for human mayus??
plots <- c("Il10")

if (species == "mouse"){
  ensembl <- useDataset("mmusculus_gene_ensembl", mart=useMart("ENSEMBL_MART_ENSEMBL"))
  BM <- getBM (attributes=c("entrezgene_id", "external_gene_name"),
               mart = ensembl, verbose = T)
id <- get(plots, org.Mm.egSYMBOL2EG)
gr <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)[id]
gene_names <- function(z) {
          names <- get(z, org.Mm.egSYMBOL)
          return(names)
          }
} else {
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
  BM <- getBM (attributes=c("entrezgene_id", "external_gene_name"),
               mart = ensembl, verbose = T)
  id <- get(plots, org.Mm.egSYMBOL2EG) #for human!
  gr <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)[id] #for human!!
  gene_names <- function(z) {
    names <- get(z, org.Mm.egSYMBOL) #for human!!
    return(names)
  }
      }

agr <- resize(gr, 100000, fix = "center")
genesinrange <- mapRangesToIds(TxDb.Mmusculus.UCSC.mm10.knownGene, agr, type = "gene")
tracks <- sapply(genesinrange[[1]][[1]], 
                 function(z) {
                   track <- geneTrack(z,TxDb.Mmusculus.UCSC.mm10.knownGene)[[1]]
                   return(track)
                 })

#Import results files from the different analysis to be included and for each replicate:
#List all the files in you data directory and import them depending on their format
file_list_bed <- list.files(path="C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/trackviewer/data",
                        pattern= "*.bed")
file_list_txt <- list.files(path="C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/trackviewer/data",
                            pattern= "*.txt")
file_list_bw <- list.files(path="C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/trackviewer/data",
                            pattern= "*.bw")

for (file in c(file_list_bed, file_list_txt)){
  x <- importScore(paste0("data/", file), 
                   format = "BED",
                   ignore.strand = T,
                   ranges = agr)
  assign(paste(strsplit(file, ".", fixed=T)[1][[1]][1]), x)
}

for (file in c(file_list_bw)){
  x <- importScore(paste0("data/", file), 
                   format = "BigWig",
                   ignore.strand = T,
                   ranges = agr)
  assign(paste(strsplit(file, ".", fixed=T)[1][[1]][1]), x)
  }


                 
#Creat functions to be used:
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
average <- function(A, B) {(A+B)/2}

  # Calculate average ATAC peaks for Tet+ samples

  average.ATAC.Tet <- Tet1_ATAC
  average.ATAC.Tet$dat2 <- average.ATAC.Tet$dat
  average.ATAC.Tet$dat <- Tet3_ATAC$dat
  average.ATAC.Tet$dat <- GRoperator(average.ATAC.Tet$dat, average.ATAC.Tet$dat2, 
                                     col="score", operator=average)
  average.ATAC.Tet$dat2 <- GRanges()
  
  # Calculate average for RNA peaks for Tet+ samples
  average.RNA.Tet <- Tet1_RNA
  average.RNA.Tet$dat <- average.RNA.Tet$dat
  average.RNA.Tet$dat2 <- Tet2_RNA$dat
  average.RNA.Tet$dat <- GRoperator(average.RNA.Tet$dat, average.RNA.Tet$dat2, 
                                    col="score", operator=average)
  average.RNA.Tet$dat2 <- GRanges()
  
  average.RNA.Tet2 <- Tet3_RNA
  average.RNA.Tet2$dat <- average_RNA.Tet2$dat
  average.RNA.Tet2$dat2 <- Tet4.RNA$dat
  average.RNA.Tet2$dat <- GRoperator(average.RNA.Tet2$dat, average.RNA.Tet2$dat2, 
                                     col="score", operator=average)
  average.RNA.Tet2$dat2 <- GRanges()
  
  av.RNA.Tet <- average.RNA.Tet
  av.RNA.Tet$dat <- av.RNA.Tet$dat
  av.RNA.Tet$dat2 <- average.RNA.Tet2$dat
  av.RNA.Tet$dat <- GRoperator(av.RNA.Tet$dat, av.RNA.Tet$dat2, 
                               col="score", operator=average)
  av.RNA.Tet$dat2 <- GRanges()
  
  # Calculate average ATAC peaks for Tconv samples
  average.ATAC.Tconv <- Tcon1.ATAC
  average.ATAC.Tconv$dat2 <- average.ATAC.Tconv$dat
  average.ATAC.Tconv$dat <- Tcon3.ATAC$dat
  average.ATAC.Tconv$dat <- GRoperator(average.ATAC.Tconv$dat, average.ATAC.Tconv$dat2, 
                                       col="score", operator=average)
  average.ATAC.Tconv$dat2 <- GRanges()
  
  # Calculate average for RNA peaks for Tconv samples
  average.RNA.Tconv <- Tconv1_RNA
  average.RNA.Tconv$dat <- average.RNA.Tconv$dat
  average.RNA.Tconv$dat2 <- Tconv2_RNA$dat
  average.RNA.Tconv$dat <- GRoperator(average.RNA.Tconv$dat, average.RNA.Tconv$dat2, 
                                      col="score", operator=average)
  average.RNA.Tconv$dat2 <- GRanges()
  
  average.RNA.Tconv2 <- Tcon3_RNA
  average.RNA.Tconv2$dat <- average.RNA.Tconv2$dat
  average.RNA.Tconv2$dat2 <- Tcon4_RNA$dat
  average.RNA.Tconv2$dat <- GRoperator(average.RNA.Tconv2$dat, average.RNA.Tconv2$dat2, 
                                       col="score", operator=average)
  average.RNA.Tconv2$dat2 <- GRanges()
  
  av.RNA.Tconv <- average.RNA.Tconv
  av.RNA.Tconv$dat <- av.RNA.Tconv$dat
  av.RNA.Tconv$dat2 <- average.RNA.Tconv2$dat
  av.RNA.Tconv$dat <- GRoperator(av.RNA.Tconv$dat, av.RNA.Tconv$dat2, 
                                 col="score", operator=average)
  av.RNA.Tconv$dat2 <- GRanges()
  
  #Create tracks for each group:
  for (i in c(1:length(pop))){
  Pattern_list<-do.call("list",mget(grep(pop[i],names(.GlobalEnv),value=TRUE)))
    o <- optimizeStyle(trackList(Pattern_list, tracks))
    assign(paste0("optSty", i), o)
    t <- o$tracks
    assign(paste0("trackList", i), t)
    v <- o$style
  assign(paste0("viewerStyle", i), v)
  }
  


  # Adjust X axis to show scale instead of chromosome ruler:
  for (i in do.call("list",mget(grep("viewerStyle",names(.GlobalEnv),value=TRUE)))){
    # Adjust X axis to show scale instead of chromosome ruler:
    setTrackViewerStyleParam(i, "xaxis", TRUE)
    setTrackViewerStyleParam(i, "margin", c(.01, .05, .01, .01))
    # Move scale in y axis to the right:
    setTrackViewerStyleParam(i, "margin", c(.08, .15, .01, .1)) # same command as before??
  }
  
  for (u in do.call("list",mget(grep("trackList",names(.GlobalEnv),value=TRUE)))){
    setTrackXscaleParam(u[[1]], "draw", FALSE)
    setTrackXscaleParam(u[[5]], "draw", TRUE) #REDUNDANT??
    setTrackXscaleParam(u[[5]], "gp", list(cex=.7))
    # Move scale in y axis to the right:
      for(i in 1:5){
      setTrackYaxisParam(u[[i]], "main", FALSE)
      }
      for(i in c(1,3:5)){
        setTrackYaxisParam(u[[i]], "gp", list(cex=.9))
      }
      # Eliminate scale for Active enhancers track:
      setTrackYaxisParam(u[[2]], "draw", FALSE)
      ## Adjust the limit of y-axis for RNA, ATAC and ChIP tracks:
      setTrackStyleParam(u[[1]], "ylim", c(0, 1))
      setTrackStyleParam(u[[2]], "ylim", c(0, 1))
      setTrackStyleParam(u[[3]], "ylim", c(0, roundUpNice(max(trackList.Tet[[3]]$dat$score,
                                                                          trackList.Tconv[[3]]$dat$score))))
      setTrackStyleParam(u[[4]], "ylim", c(0, roundUpNice(max(trackList.Tet[[4]]$dat$score,
                                                                          trackList.Tconv[[4]]$dat$score))))
      setTrackStyleParam(u[[5]], "ylim", c(0, roundUpNice(max(trackList.Tet[[5]]$dat$score,
                                                                          trackList.Tconv[[5]]$dat$score))))
      for(i in 1:5){
        setTrackStyleParam(u[[i]], "marginBottom", .1)
      }
      # For each transcript, the transcript name can be put on the upstream or downstream 
      #of the transcript
      for(i in 6:length(u)){
        setTrackStyleParam(u[[i]], "ylabpos", "upstream")
        setTrackStyleParam(u[[i]], "ylabgp", list(cex=1))
      }
      # Adjust the label of y-axis (colour and size)
      for (i in 1:5){
                setTrackStyleParam(u[[i]], "ylabgp", list(cex=.8, col="black")) 
      }
      # Adjust track colors:
      setTrackStyleParam(u[[1]], "color", "darkgreen")
      setTrackStyleParam(u[[2]], "color", "red")
      setTrackStyleParam(u[[3]], "color", "red")
      setTrackStyleParam(u[[4]], "color", "purple3")
      setTrackStyleParam(u[[5]], "color", "darkblue")
      for(i in 6:length(u)){
        setTrackStyleParam(u[[i]], "color", "grey35")
      }
      # Change the track names
      names(tracks) <- sapply(names(tracks),gene_names) #check function earlier
                              
      names(u) <- c("DIFF METH", "DIFF ACT ENH", "H3K27ac", 
                                "ATAC", "RNA", names(tracks))
      # Adjust the track height
      setTrackStyleParam(u[[2]], "height", 0.04)
      for(i in c(1,3:5)){
        setTrackStyleParam(u[[i]], "height", 0.170)
      }
      for(i in 6:length(u)){
        setTrackStyleParam(u[[i]], "height", 0.04)
      }
      
    }
  
 
 ##script below, general or specific for IL-10..? 
  
  Il10.enh.DMR <- Il10.enh.DMR[order(as.numeric(gsub("chr", "", Il10.enh.DMR$Chr)), 
                                    as.numeric(Il10.enh.DMR$Start),
                                    decreasing = F, na.last = T), ]
  
  vp <- viewTracks(trackList.Tet, gr=agr, viewerStyle=viewerStyle.Tet)
  addGuideLine(c(start(gr), end(gr)), vp=vp)
#  addGuideLine(c(Il10.enh.DMR$Start[1], Il10.enh.DMR$End[1]), vp=vp)
#  addGuideLine(c(Il10.enh.DMR$Start[2], Il10.enh.DMR$End[2]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[3], Il10.enh.DMR$End[3]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[4], Il10.enh.DMR$End[4]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[5], Il10.enh.DMR$End[5]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[6], Il10.enh.DMR$End[6]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[7], Il10.enh.DMR$End[7]), vp=vp)
  addGuideLine(c(Il10.enh.DMR$Start[8], Il10.enh.DMR$End[8]), vp=vp)

methy <- read.table("/Users/patri/Documents/LAB/TESIS DOCTORAL 2015-2020/TR1 PROJECT/2019_05_Methylome_Tet_Tconv/Methylome results_filtered/DMR_Tet_Tconv/CG.group.site.methy_ratio.txt",
                    sep = "\t", dec = ".", header = T, quote = "")

Tet <- GRanges(seqnames = methy$chr, 
               ranges = IRanges(start = methy$pos,
                                width = 1), 
               strand = NULL,
               #We cross here columns because we assume groups were crossed
               #in the analysis!!
               score = as.integer(methy$Tet*100))

Tconv <- GRanges(seqnames = methy$chr, 
                 ranges = IRanges(start = methy$pos,
                                  width = 1), 
                 strand = NULL,
                 #We cross here columns because we assume groups were crossed
                 #in the analysis!!
                 score = as.integer(methy$Tconv*100))

#Lolliplots overlapping Tet and Tconv:
pdf("Differential single Cs in Il10 associated active enhancers.pdf", 
    width = 12, height = 6)
vp <- viewTracks(trackList.Tet, gr=agr, viewerStyle=viewerStyle.Tet)
addGuideLine(c(start(gr), end(gr)), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[1], Il10.enh.DMR$End[1]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[2], Il10.enh.DMR$End[2]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[3], Il10.enh.DMR$End[3]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[4], Il10.enh.DMR$End[4]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[5], Il10.enh.DMR$End[5]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[6], Il10.enh.DMR$End[6]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[7], Il10.enh.DMR$End[7]), vp=vp)
addGuideLine(c(Il10.enh.DMR$Start[8], Il10.enh.DMR$End[8]), vp=vp)

for (i in 1:nrow(Il10.enh.DMR)) {
grl <- c(Tet, Tconv)
ex.ac.enh <- GRanges(seqnames = Il10.enh.DMR$Chr[i], 
                     ranges = paste0(Il10.enh.DMR$Start[i], "-", Il10.enh.DMR$End[i]), 
                     strand = NULL)
  gr <- ex.ac.enh
  Tet$color <- "#FFD700"
  Tet$border <- "gray30"
  Tconv$color <- "#DB7575"
  Tconv$border <- "gray30"
  xaxis <- c(seq(from = start(gr), 
                 to = end(gr), 
                 by = 500))
  yaxis <- c(0, 100)
  legends <- list(list(labels=c("Tet", "Tconv"), 
                       fill=c("#FFD700", "#DB7575"),
                       col = c("gray30","gray30")))
  
  lolliplot(grl, gr, type = "circle", 
            xaxis= xaxis,
            yaxis= yaxis,
            legend= legends,
            ylab= paste0(seqnames(gr), ":", start(gr), "-", end(gr)))
}
dev.off()
  
