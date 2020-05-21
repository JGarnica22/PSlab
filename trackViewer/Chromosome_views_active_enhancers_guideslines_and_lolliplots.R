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
library(biomaRt)

# Set your working directory (the project you are working in):
setwd("/home/jgarnica/R/trackviewer")

# Indicate populations to work with
pop <- c("Tconv", "Tet")

# Indicate species working with (mouse or human)
species <- "mouse"

# Indicate the genes to visualize
plots <- c("Il10", "Il21", "Mapkapk2", "Il19" )
# or read a txt file containing a list with all the plots

# Create functions to be used:
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
average <- function(A, B) {(A+B)/2}
average3 <- function (A,B) {(A+B)/3}


if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
}
BM <- getBM (attributes=c("entrezgene_id", "external_gene_name"),
             mart = ensembl, verbose = T)

#Load methylation comparison of samples
CG.methy <- read.table(paste0("data/", list.files(path="/home/jgarnica/R/trackviewer/data", pattern= "*CG.*")), 
                       sep = "\t", dec = ".", header = T, quote = "")

for (xy in plots) {
  id <- get(xy, org.SYMBOL2EG)
  gr <- genes(TxDb)[id]
  agr <- if (width(gr)<120000) {
    resize(gr, 100000, fix = "center")
  } else {
    resize(gr, width(gr)+20000, fix = "center")
  }
  # agr.gene <- if (width(gr)<100000) {
  #   resize(gr, 300000, fix = "center")
  # } else {
  #   resize(gr, width(gr)+50000, fix = "center")
  # }
  genesinrange <- mapRangesToIds(TxDb, agr, type = "gene", ignore.strand = T)
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
  file_list_DMR <- list.files(path="/home/jgarnica/R/trackviewer/data",
                              pattern= "*DMR*")
  
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
  for (l in 1:length(pop)){
    for (s in c("ATAC", "RNA", "ChIP", "enh", "methy")){  
      for (i in 1:length(grep(paste0(pop[l], ".*" , s), names(.GlobalEnv),value=TRUE))){
        me <- eval(as.symbol(grep(paste0(pop[l], ".*" , s), names(.GlobalEnv),value=TRUE)[i]))
        assign(paste0(s, ".", comp[l], i), me)
      }
    }
    
    #indicate experiments with multiple replicates
    for (n in c("ATAC", "RNA")){
      #assign all objects to be used for operations
      for(w in 1:length(grep(paste0(n, "*.",comp[l], "."), names(.GlobalEnv),value=TRUE))){
        a <- eval(as.symbol(grep(paste0(n, ".",comp[l], "."), names(.GlobalEnv),value=TRUE)[w]))
        assign(paste0("mit",w), a)
      }
      
      if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 2) {
        avr <- mit1
        avr$dat <- avr$dat
        avr$dat2 <- mit2$dat
        avr$dat <- GRoperator(avr$dat, avr$dat2, 
                              col="score", operator=average)
        avr$dat2 <- GRanges()
        assign(paste0("arch", n, ".", comp[l]), avr)
      }
      
      if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 3){
        avr <- mit1
        avr$dat <- avr$dat
        avr$dat2 <- mit2$dat
        avr$dat <- GRoperator(avr$dat, avr$dat2, 
                              col="score", operator="+")
        avr$dat2 <- mit3$dat
        avr$dat <- GRoperator(avr$dat, avr$dat2, 
                              col="score", operator=average3)
        avr$dat2 <- GRanges()
        assign(paste0("arch", n, ".", comp[l]), avr)
      }
      
      if (length(grep("mit", names(.GlobalEnv),value=TRUE)) == 4){
        avrg <- mit1
        avrg$dat <- avrg$dat
        avrg$dat2 <- mit2$dat
        avrg$dat <- GRoperator(avrg$dat, avrg$dat2, 
                               col="score", operator=average)
        avrg$dat2 <- GRanges()
        
        avrg2 <- mit3
        avrg2$dat <- avrg2$dat
        avrg2$dat2 <- mit4$dat
        avrg2$dat <- GRoperator(avrg2$dat, avrg2$dat2, 
                                col="score", operator=average)
        avrg2$dat2 <- GRanges()
        
        avr <- avrg
        avr$dat <- avr$dat
        avr$dat2 <- avrg2$dat
        avr$dat <- GRoperator(avr$dat, avr$dat2, 
                              col="score", operator=average)
        avr$dat2 <- GRanges()
        assign(paste0("arch", n, ".", comp[l]), avr)
      }
      rm(list=grep("mit", names(.GlobalEnv), value = T))
    }
  }
  
  
  # Create tracks for each group, add "arch" for replicates:
  show.tracks <- c(mget(grep("methy.",names(.GlobalEnv),value=TRUE, fixed = T)),
                   mget(grep("enh.",names(.GlobalEnv),value=TRUE, fixed = T)),
                   mget(grep("ChIP.",names(.GlobalEnv),value=TRUE, fixed = T)),
                   mget(grep("archATAC.",names(.GlobalEnv),value=TRUE, fixed = T)),
                   mget(grep("archRNA.",names(.GlobalEnv),value=TRUE, fixed = T)))
  
  # Change the track names
  names(tracks) <- sapply(names(tracks), 
                          function(z) {
                            names <- get(z, org.Mm.egSYMBOL)
                            return(names)
                          })
  #Prepare Granges objects for lolliplots over methylation status on Active enhancer mapping
  for (gh in 1:length(pop)){
    y <- GRanges(seqnames = CG.methy$chr, ranges = IRanges(start = CG.methy$pos,
                                                           width = 1), strand = NULL,
                 #We cross here columns because we assume groups were crossed
                 #in the analysis!!
                 score = as.integer(CG.methy[,gh+2]*100))
    assign(paste0("lolli.", pop[gh]), y)
  }
  grl <- c(mget(grep("lolli.*",names(.GlobalEnv),value=TRUE)))
  colors <- c("#FFD700", "#DB7575", "#81C784")[1:length(pop)]
  grlo <- GRanges()
  for (p in 1:length(grl)){
    grl[[p]]$color <- colors[p]
    grl[[p]]$border <- "gray30"
    grlo <- c(grlo, grl[[p]])
  }
  
  pdf(file = paste0("figs/",xy, "_active_enhancers.pdf"), width = 10, height = 6)
  for (m in 1:length(pop)){
    Pattern_list<-do.call("list",mget(grep(comp[m],names(show.tracks),value=TRUE)))
    o <- optimizeStyle(trackList(Pattern_list, tracks))
    assign(paste0("optSty", m), o)
    u <- o$tracks
    assign(paste0("trackList", m), t)
    t <- o$style
    assign(paste0("viewerStyle", m), u)
    
    ##Generate active enhancers guidelines
    #Read table of previously elaborated active enhancers map and prepare dataframe
    hun100kb <- read.table(paste0("data/", grep(pop[m], file_list_DMR, value = T)))
    hun100kb <- hun100kb[, c(1:3, 9, 7)]
    names(hun100kb) <- c("Chr", "Start", "End", "Strand", "EntrezID")
    hun100kb <- merge(hun100kb, BM, by.x = "EntrezID", by.y = "entrezgene_id")
    hun100kb <- hun100kb[, c(-1)]
    names(hun100kb) <- c("Chr", "Start", "End", "Strand", "Gene")
    
    #select only active enhancers for the gene to be represented
    xy.enh.DMR <-  hun100kb[which(hun100kb$Gene==xy),]
    xy.enh.DMR <- xy.enh.DMR[order(as.numeric(gsub("chr", "", xy.enh.DMR$Chr)), 
                                   as.numeric(xy.enh.DMR$Start),
                                   decreasing = F, na.last = T), ]
    
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
      grx <- genes(TxDb)[as.data.frame(genesinrange)[i-5,1]]
      if (start(ranges(grx)) <= start(ranges(agr))){
          setTrackStyleParam(u[[i]], "ylabpos", "left")
      } 
      if(end(ranges(grx)) >= end(ranges(agr))){
        setTrackStyleParam(u[[i]], "ylabpos", "right")
      } else {
        setTrackStyleParam(u[[i]], "ylabpos", "upstream")
      } 
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
    } else {
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
    #add guidelines to mark inferred active enhancers
    if (nrow(xy.enh.DMR)!=0){
      for (c in 1:nrow(xy.enh.DMR)){
        addGuideLine(c(xy.enh.DMR$Start[c], xy.enh.DMR$End[c]), vp=vp)
       }
    }
  }
    dev.off()
    
    for (mi in 1:length(pop)){
    #Read table of previously elaborated active enhancers map and prepare dataframe
    hun100kb <- read.table(paste0("data/", grep(pop[mi], file_list_DMR, value = T)))
    hun100kb <- hun100kb[, c(1:3, 9, 7)]
    names(hun100kb) <- c("Chr", "Start", "End", "Strand", "EntrezID")
    hun100kb <- merge(hun100kb, BM, by.x = "EntrezID", by.y = "entrezgene_id")
    hun100kb <- hun100kb[, c(-1)]
    names(hun100kb) <- c("Chr", "Start", "End", "Strand", "Gene")
    
    #select only active enhancers for the gene to be represented
    xy.enh.DMR <-  hun100kb[which(hun100kb$Gene==xy),]
    xy.enh.DMR <- xy.enh.DMR[order(as.numeric(gsub("chr", "", xy.enh.DMR$Chr)), 
                                   as.numeric(xy.enh.DMR$Start),
                                   decreasing = F, na.last = T), ]
    if (nrow(xy.enh.DMR)!=0){
      pdf(paste0("figs/Lolliplot_",xy, "_", paste0(pop, collapse ="_"), "_on_", pop[mi], "_AE", ".pdf"), width = 12, height = 6)
      for (ic in 1:nrow(xy.enh.DMR)){
                #Do also the lolliplots for nucleotid-specific metilation status on the active enhancers
        grpi <- GRanges(seqnames = xy.enh.DMR$Chr[ic], 
                        ranges = paste0(xy.enh.DMR$Start[ic], "-", xy.enh.DMR$End[ic]), 
                        strand = NULL)
        xaxis <- c(seq(from = start(gr), 
                       to = end(gr), 
                       by = 500))
        yaxis <- c(0, 100)
        legends <- list(list(labels=c(pop[2],pop[1]), 
                             fill= colors,
                             col = c("gray30","gray30")))
        lolliplot(grlo, grpi, type = "circle", 
                  xaxis= xaxis,
                  yaxis= yaxis,
                  legend= legends,
                  ylab= paste0(xy, "-", seqnames(grpi), ":", start(grpi), "-", end(grpi)))
      }
      dev.off()
      }
    }
  
}
