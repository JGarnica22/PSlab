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
setwd("/home/jgarnica/R/trackviewer")

# Indicate populations to work with
pop <- c("Tconv", "Tet")

# Indicate species working with (mouse or human)
species <- "mouse"

# Indicate the genes to visualize
plots <- c("Il10", "Il21")
# or read a txt file containing a list with all the plots

# Create functions to be used:
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
average <- function(A, B) {(A+B)/2}
average3 <- function (A,B) {(A+B*2)/3}


if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
}

for (xy in plots) {
  id <- get(xy, org.SYMBOL2EG)
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
  #IMPORTANT: name files always as type_of_sample(Tet, Tconv, TFH...)+number of replicate_tecnique (ATAC, RNA...), e.g. Tet3_ATAC, Tet1_RNA...
  file_list_bed <- list.files(path="/home/jgarnica/R/trackviewer/data",
                              pattern= "*.bed")
  file_list_bw <- list.files(path="/home/jgarnica/R/trackviewer/data",
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
  
  
  # Calculate average peaks for all assays and samples
  for (l in c(1:length(pop))){
    for (s in c("ATAC", "RNA", "ChIP", "enh", "methy")){  
      for (i in 1:length(grep(paste0(pop[l], ".*" , s), names(.GlobalEnv),value=TRUE))){
        me <- eval(as.symbol(grep(paste0(pop[l], ".*" , s), names(.GlobalEnv),value=TRUE)[i]))
        assign(paste0(s, ".", comp[l], i), me)
      }
    }
    #indicate experiments with multiple replicates
    for (n in c("ATAC", "RNA")){
      #assign all objects to be used for operations
      for(w in 1:length(grep(paste0(n, "*.",comp[l]), names(.GlobalEnv),value=TRUE))){
        a <- eval(as.symbol(grep(paste0(n, "*.",comp[l]), names(.GlobalEnv),value=TRUE)[w]))
        assign(paste0("mit",w), a)
      }
      #Operate them depending on the number of ojbects, final objects will be names as "arch" and intermediate "avg"
        for (j in seq(1, length(grep("mit", names(.GlobalEnv),value=TRUE)), by=2)) {
          for (g in seq(2, length(grep("mit", names(.GlobalEnv),value=TRUE)), by=2)){
            
            if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 2) {
              avr <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[j]))
              avr$dat2 <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[j]))$dat
              avr$dat <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[g]))$dat
              avr$dat <- GRoperator(avr$dat, avr$dat2, 
                                    col="score", operator=average)
              avr$dat2 <- GRanges()
              assign(paste0("arch" ,n, ".", comp[l]), avr)
            } 
            if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 4){
              avr <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[j]))
              avr$dat2 <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[j]))$dat
              avr$dat <- eval(as.symbol(grep("mit", names(.GlobalEnv),value=TRUE)[g]))$dat
              avr$dat <- GRoperator(avr$dat, avr$dat2, 
                                    col="score", operator=average)
              avr$dat2 <- GRanges()
              assign(paste0("avg" ,j), avr)
              }
              
              if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 3){
                avr <- mit1
                avr$dat2 <- mit1$dat
                avr$dat <- mit2$dat
                avr$dat <- GRoperator(avr$dat, avr$dat2, 
                                      col="score", operator=average)
                avr$dat2 <- GRanges()
                assign(paste0("avgtr" ,n, ".", comp[l]), avr)
                #do average using average3 functions for triplicates
                  avrg <- eval(as.symbol(grep("avgtr*", names(.GlobalEnv),value=TRUE)[j]))
                  avrg$dat2 <- mit3$dat
                  avrg$dat <- eval(as.symbol(grep(paste0("avgtr"), names(.GlobalEnv),value=TRUE)))$dat
                  avrg$dat <- GRoperator(avrg$dat, avrg$dat2, 
                                        col="score", operator=average3)
                  avrg$dat2 <- GRanges()
                  assign(paste0("arch" ,n, ".", comp[l]), avrg)
                    } 
                  }
                }
              #average of averages for 4 samples, need to be outside j - g loop
            if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 4){
                avrg <- avg1
                avrg$dat2 <- avg1$dat
                avrg$dat <- avg3$dat
                avrg$dat <- GRoperator(avrg$dat, avrg$dat2, 
                                      col="score", operator=average)
                avrg$dat2 <- GRanges()
                assign(paste0("arch" ,n, ".", comp[l]), avrg)
              }
          }
        
      }
  
 #former way
  av.ATAC.Ctl <- ATAC.Ctl1
  av.ATAC.Ctl$dat2 <- av.ATAC.Ctl$dat
  av.ATAC.Ctl$dat <- ATAC.Ctl2$dat
  av.ATAC.Ctl$dat <- GRoperator(av.ATAC.Ctl$dat, av.ATAC.Ctl$dat2, 
                                col="score", operator=average)
  av.ATAC.Ctl$dat2 <- GRanges()

  archATAC.Ctl

  
  average.RNA.Ctl <- RNA.Ctl1
  average.RNA.Ctl$dat <- average.RNA.Ctl$dat
  average.RNA.Ctl$dat2 <- RNA.Ctl2$dat
  average.RNA.Ctl$dat <- GRoperator(average.RNA.Ctl$dat, average.RNA.Ctl$dat2, 
                                    col="score", operator=average)
  average.RNA.Ctl$dat2 <- GRanges()
  avg1
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
  archRNA.Ctl
  
  av.ATAC.Smp <- ATAC.Smp1
  av.ATAC.Smp$dat2 <- av.ATAC.Smp$dat
  av.ATAC.Smp$dat <- ATAC.Smp2$dat
  
  av.ATAC.Smp$dat <- GRoperator(av.ATAC.Smp$dat, av.ATAC.Smp$dat2, 
                                col="score", operator=average)
  av.ATAC.Smp$dat2 <- GRanges()
  archATAC.Smp
  
  average.RNA.Smp <- RNA.Smp1
  average.RNA.Smp$dat <- average.RNA.Smp$dat
  average.RNA.Smp$dat2 <- RNA.Smp2$dat
  average.RNA.Smp$dat <- GRoperator(average.RNA.Smp$dat, average.RNA.Smp$dat2, 
                                    col="score", operator=average)
  average.RNA.Smp$dat2 <- GRanges()
  avg1
  
  average.RNA.Smp2 <- RNA.Smp3
  average.RNA.Smp2$dat <- average.RNA.Smp2$dat
  average.RNA.Smp2$dat2 <- RNA.Smp4$dat
  average.RNA.Smp2$dat <- GRoperator(average.RNA.Smp2$dat, average.RNA.Smp2$dat2, 
                                     col="score", operator=average)
  average.RNA.Smp2$dat2 <- GRanges()
  avg3
  av.RNA.Smp <- average.RNA.Smp
  av.RNA.Smp$dat <- av.RNA.Smp$dat
  av.RNA.Smp$dat2 <- average.RNA.Smp2$dat
  av.RNA.Smp$dat <- GRoperator(av.RNA.Smp$dat, av.RNA.Smp$dat2, 
                               col="score", operator=average)
  av.RNA.Smp$dat2 <- GRanges()
  archRNA.Smp
  
  
  # Create tracks for each group, add "arch" for replicates:
  show.tracks <- c(mget(grep("methy.",names(.GlobalEnv),value=TRUE)),
                   mget(grep("enh.",names(.GlobalEnv),value=TRUE)),
                   mget(grep("ChIP.",names(.GlobalEnv),value=TRUE)),
                   mget(grep("archATAC.",names(.GlobalEnv),value=TRUE)),
                   mget(grep("archRNA.",names(.GlobalEnv),value=TRUE)))
  
  # Change the track names
  names(tracks) <- sapply(names(tracks), 
                          function(z) {
                            names <- get(z, org.Mm.egSYMBOL)
                            return(names)
                          })
  
  pdf(file = paste0("figs/",xy, ".pdf"), width = 10, height = 6)
  for (m in c(1:length(pop))){
    Pattern_list<-do.call("list",mget(grep(comp[m],names(show.tracks),value=TRUE)))
    o <- optimizeStyle(trackList(Pattern_list, tracks))
    assign(paste0("optSty", m), o)
    u <- o$tracks
    assign(paste0("trackList", m), t)
    t <- o$style
    assign(paste0("viewerStyle", m), u)
    
    # Adjust X axis to show scale instead of chromosome ruler:
    setTrackViewerStyleParam(t, "xaxis", FALSE)
    setTrackViewerStyleParam(t, "margin", c(.01, .05, .01, .01))
    
    setTrackXscaleParam(u[[1]], "draw", FALSE)
    setTrackXscaleParam(u[[5]], "draw", TRUE) #REDUNDANT??
    setTrackXscaleParam(u[[5]], "gp", list(cex=.7))
    # Move scale in y axis to the right:
    setTrackViewerStyleParam(t, "margin", c(.08, .15, .01, .1))
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
    setTrackStyleParam(u[[3]], "ylim", c(0, roundUpNice(max(optimizeStyle(trackList(do.call("list",mget(grep(comp[1],names(show.tracks),value=TRUE))), tracks))$tracks[[3]]$dat$score,
                                                            optimizeStyle(trackList(do.call("list",mget(grep(comp[2],names(show.tracks),value=TRUE))), tracks))$tracks[[3]]$dat$score))))
    setTrackStyleParam(u[[4]], "ylim", c(0, roundUpNice(max(optimizeStyle(trackList(do.call("list",mget(grep(comp[1],names(show.tracks),value=TRUE))), tracks))$tracks[[4]]$dat$score,
                                                            optimizeStyle(trackList(do.call("list",mget(grep(comp[2],names(show.tracks),value=TRUE))), tracks))$tracks[[4]]$dat$score))))
    setTrackStyleParam(u[[5]], "ylim", c(0, roundUpNice(max(optimizeStyle(trackList(do.call("list",mget(grep(comp[1],names(show.tracks),value=TRUE))), tracks))$tracks[[5]]$dat$score,
                                                            optimizeStyle(trackList(do.call("list",mget(grep(comp[2],names(show.tracks),value=TRUE))), tracks))$tracks[[5]]$dat$score))))
    
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
    
    names(u) <- c("DIFF METH", "DIFF ACT ENH", "H3K27ac", 
                  "ATAC", "RNA", names(tracks))
    
    # Adjust the track height, add conditionals in case many genes are added
    if (length(u) <= 10){
    setTrackStyleParam(u[[2]], "height", 0.04)
    for(i in c(1,3:5)){
      setTrackStyleParam(u[[i]], "height", 0.170)
    }
    for(i in 6:length(u)){
      setTrackStyleParam(u[[i]], "height", 0.04)
    }
    } if (length(u) > 10){
      setTrackStyleParam(u[[2]], "height", 0.04)
      for(i in c(1,3:5)){
        setTrackStyleParam(u[[i]], "height", 0.170-((length(u)-10)*0.01))
      }
      for(i in 6:length(u)){
        setTrackStyleParam(u[[i]], "height", 0.04)
      } 
      
    } 
    
  
    if (m == 1){
      vp <- viewTracks(u, gr=agr, viewerStyle=t, newpage = F)
    } else {
      vp <- viewTracks(u, gr=agr, viewerStyle=t, newpage = T)
    }
    addGuideLine(c(start(gr), end(gr)), vp=vp)
    grid.text(pop[m], x=0.2, y=.9, just="top", 
              gp=gpar(cex=1.5, fontface="bold"))
    
  }
  dev.off()
}
