#This script is to visualize genomic data by trackViewer
#It was created for R 3.6.3 version (2020-05-03)
#Copyright (C) 2020  Patricia Solé Sánchez and Josep Garnica Caparrós
#####################################################################

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

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/Functionals/trackviewer")

# Indicate populations to work with
pop <- c("Tconv", "Tet")

# Indicate species working with (mouse or human)
species <- "mouse"

# Indicate the genes to visualize
# plots <- c("Il10")
# or read a txt file containing a list with all the plots

# Create functions to be used:
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
average <- function(A, B) {(A+B)/2}


if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
}


plot <- "Il10"

#for (plot in plots) {
id <- get(plot, org.SYMBOL2EG)
gr <- genes(TxDb)[id]
agr <- if (width(gr)<100000) {
  resize(gr, 100000, fix = "center")
} else {
  resize(gr, width(gr)+20000, fix = "center")
}
genesinrange <- mapRangesToIds(TxDb, agr, type = "gene")
tracks <- sapply(genesinrange[[1]][[1]], 
                 function(z) {
                   track <- geneTrack(z,TxDb)[[1]]
                   return(track)
                 })

#Import results files from the different analysis to be included and for each replicate:
#List all the files in you data directory and import them depending on their format
file_list_bed <- list.files(path="/Users/patri/Documents/LAB/TESIS DOCTORAL 2015-2020/INTERNSHIP/8_trackViewer/data",
                            pattern= "*.bed")
file_list_bw <- list.files(path="/Users/patri/Documents/LAB/TESIS DOCTORAL 2015-2020/INTERNSHIP/8_trackViewer/data",
                           pattern= "*.bw")

for (file in c(file_list_bed)){
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

comp <- c("Ctl", "Smp")


# Calculate average ATAC peaks for pop[1] samples
for (i in 1:length(grep(paste0(pop[1],".*ATAC"),names(.GlobalEnv),value=TRUE))){
  ATAC.Ctl1 <- eval(as.symbol(grep(paste0(pop[1],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i]))
  ATAC.Ctl2 <- eval(as.symbol(grep(paste0(pop[1],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+1]))
  ATAC.Ctl3 <- eval(as.symbol(grep(paste0(pop[1],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+2]))
  ATAC.Ctl4 <- eval(as.symbol(grep(paste0(pop[1],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+3]))
}

av.ATAC.Ctl <- ATAC.Ctl1
av.ATAC.Ctl$dat2 <- av.ATAC.Ctl$dat
av.ATAC.Ctl$dat <- ATAC.Ctl2$dat

av.ATAC.Ctl$dat <- GRoperator(av.ATAC.Ctl$dat, av.ATAC.Ctl$dat2, 
                              col="score", operator=average)
av.ATAC.Ctl$dat2 <- GRanges()

# Calculate average for RNA peaks for pop[1] samples
for (i in 1:length(grep(paste0(pop[1],".*RNA"),names(.GlobalEnv),value=TRUE))){
  RNA.Ctl1 <- eval(as.symbol(grep(paste0(pop[1],".*RNA"),names(.GlobalEnv),,value=TRUE)[i]))
  RNA.Ctl2 <- eval(as.symbol(grep(paste0(pop[1],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+1]))
  RNA.Ctl3 <- eval(as.symbol(grep(paste0(pop[1],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+2]))
  RNA.Ctl4 <- eval(as.symbol(grep(paste0(pop[1],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+3]))
}

average.RNA.Ctl <- RNA.Ctl1
average.RNA.Ctl$dat <- average.RNA.Ctl$dat
average.RNA.Ctl$dat2 <- RNA.Ctl2$dat
average.RNA.Ctl$dat <- GRoperator(average.RNA.Ctl$dat, average.RNA.Ctl$dat2, 
                                  col="score", operator=average)
average.RNA.Ctl$dat2 <- GRanges()

average.RNA.Ctl2 <- RNA.Ctl3
average.RNA.Ctl2$dat <- average.RNA.Ctl2$dat
average.RNA.Ctl2$dat2 <- RNA.Ctl4$dat
average.RNA.Ctl2$dat <- GRoperator(average.RNA.Ctl2$dat, average.RNA.Ctl2$dat2, 
                                   col="score", operator=average)
average.RNA.Ctl2$dat2 <- GRanges()

av.RNA.Ctl <- average.RNA.Ctl
av.RNA.Ctl$dat <- av.RNA.Ctl$dat
av.RNA.Ctl$dat2 <- average.RNA.Ctl2$dat
av.RNA.Ctl$dat <- GRoperator(av.RNA.Ctl$dat, av.RNA.Ctl$dat2, 
                             col="score", operator=average)
av.RNA.Ctl$dat2 <- GRanges()


# Calculate average ATAC peaks for pop[2] samples
for (i in 1:length(grep(paste0(pop[2],".*ATAC"),names(.GlobalEnv),value=TRUE))){
  ATAC.Smp1 <- eval(as.symbol(grep(paste0(pop[2],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i]))
  ATAC.Smp2 <- eval(as.symbol(grep(paste0(pop[2],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+1]))
  ATAC.Smp3 <- eval(as.symbol(grep(paste0(pop[2],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+2]))
  ATAC.Smp4 <- eval(as.symbol(grep(paste0(pop[2],".*ATAC"),names(.GlobalEnv),,value=TRUE)[i+3]))
}

av.ATAC.Smp <- ATAC.Smp1
av.ATAC.Smp$dat2 <- av.ATAC.Smp$dat
av.ATAC.Smp$dat <- ATAC.Smp2$dat

av.ATAC.Smp$dat <- GRoperator(av.ATAC.Smp$dat, av.ATAC.Smp$dat2, 
                              col="score", operator=average)
av.ATAC.Smp$dat2 <- GRanges()

# Calculate average for RNA peaks for pop[2] samples
for (i in 1:length(grep(paste0(pop[2],".*RNA"),names(.GlobalEnv),value=TRUE))){
  RNA.Smp1 <- eval(as.symbol(grep(paste0(pop[2],".*RNA"),names(.GlobalEnv),,value=TRUE)[i]))
  RNA.Smp2 <- eval(as.symbol(grep(paste0(pop[2],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+1]))
  RNA.Smp3 <- eval(as.symbol(grep(paste0(pop[2],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+2]))
  RNA.Smp4 <- eval(as.symbol(grep(paste0(pop[2],".*RNA"),names(.GlobalEnv),,value=TRUE)[i+3]))
}

average.RNA.Smp <- RNA.Smp1
average.RNA.Smp$dat <- average.RNA.Smp$dat
average.RNA.Smp$dat2 <- RNA.Smp2$dat
average.RNA.Smp$dat <- GRoperator(average.RNA.Smp$dat, average.RNA.Smp$dat2, 
                                  col="score", operator=average)
average.RNA.Smp$dat2 <- GRanges()

average.RNA.Smp2 <- RNA.Smp3
average.RNA.Smp2$dat <- average.RNA.Smp2$dat
average.RNA.Smp2$dat2 <- RNA.Smp4$dat
average.RNA.Smp2$dat <- GRoperator(average.RNA.Smp2$dat, average.RNA.Smp2$dat2, 
                                   col="score", operator=average)
average.RNA.Smp2$dat2 <- GRanges()

av.RNA.Smp <- average.RNA.Smp
av.RNA.Smp$dat <- av.RNA.Smp$dat
av.RNA.Smp$dat2 <- average.RNA.Smp2$dat
av.RNA.Smp$dat <- GRoperator(av.RNA.Smp$dat, av.RNA.Smp$dat2, 
                             col="score", operator=average)
av.RNA.Smp$dat2 <- GRanges()


ChIP.Ctl <- eval(as.symbol(grep(paste0(pop[1],".*ChIP"),names(.GlobalEnv),,value=TRUE)))
ChIP.Smp <- eval(as.symbol(grep(paste0(pop[2],".*ChIP"),names(.GlobalEnv),,value=TRUE)))

Enh.Ctl <- eval(as.symbol(grep(paste0(pop[1],".*enh"),names(.GlobalEnv),,value=TRUE)))
Enh.Smp <- eval(as.symbol(grep(paste0(pop[2],".*enh"),names(.GlobalEnv),,value=TRUE)))

Methy.Ctl <- eval(as.symbol(grep(paste0(pop[1],".*methy"),names(.GlobalEnv),,value=TRUE)))
Methy.Smp <- eval(as.symbol(grep(paste0(pop[2],".*methy"),names(.GlobalEnv),,value=TRUE)))


# Create tracks for each group:
show.tracks <- c(mget(grep("Methy.",names(.GlobalEnv),value=TRUE)),
                 mget(grep("Enh.",names(.GlobalEnv),value=TRUE)),
                 mget(grep("ChIP.",names(.GlobalEnv),value=TRUE)),
                 mget(grep("av.ATAC.",names(.GlobalEnv),value=TRUE)),
                 mget(grep("av.RNA.",names(.GlobalEnv),value=TRUE)))

# Change the track names
names(tracks) <- sapply(names(tracks), 
                        function(z) {
                          names <- get(z, org.Mm.egSYMBOL)
                          return(names)
                        })


for (i in c(1:length(pop))){
  Pattern_list<-do.call("list",mget(grep(comp[i],names(show.tracks),value=TRUE)))
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
  setTrackStyleParam(u[[3]], "ylim", c(0, roundUpNice(max(trackList1[[3]]$dat$score,
                                                          trackList2[[3]]$dat$score))))
  setTrackStyleParam(u[[4]], "ylim", c(0, roundUpNice(max(trackList1[[4]]$dat$score,
                                                          trackList2[[4]]$dat$score))))
  setTrackStyleParam(u[[5]], "ylim", c(0, roundUpNice(max(trackList1[[5]]$dat$score,
                                                          trackList2[[5]]$dat$score))))
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
  
  n <- grep(strsplit(names(u[1]), ".", fixed = T)[[1]][2], comp)
  
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
  vp <- viewTracks(u, gr=agr, viewerStyle=viewerStyle1)
  addGuideLine(c(start(gr), end(gr)), vp=vp)
  grid.text(pop[n],
            x=.5, y=.9, just="top", 
            gp=gpar(cex=1.5, fontface="bold"))
}
