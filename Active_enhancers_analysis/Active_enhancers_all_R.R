#This script is to infer active enhancer regions based on ATAC anc ChIP-seq data
#It was created for R 3.6.3 version (2020-05-22)
#Copyright (C) 2020  Patricia Solé Sánchez and Josep Garnica Caparrós
#####################################################################
# Check if required packages are installed, if not install:
cran.packages <- c("Cairo", "stringr", "tidyr", "snakecase", "plyr", "writexl")
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
                   "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "ggplot")
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
library(stringr)
library(ggplot2)
library(gridExtra)

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
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
}
BM <- getBM (attributes=c("entrezgene_id", "external_gene_name"),
             mart = ensembl, verbose = T)

#Load all files needed
#Load DMR file between two samples:
DMR <- read.table("data/DMR.txt",
                  sep = "\t", quote = "",
                  dec = ".", header = T, na.strings = T)
names(DMR) <- c("Chr", "Start", "End", pop[1], pop[2])
DMR <- DMR[order(as.numeric(gsub("chr", "", DMR$Chr)), 
                 as.numeric(DMR$Start),
                 decreasing = F, na.last = T), ]

#Add to your data directory the DESeq2 file comparing you samples to analyse

DESeq2 <- read.table (file = paste0("data/", list.files(path=paste0(getwd(),"/data"), pattern= "DESeq2_")),
                      sep = "\t", quote = "", dec = ".", header=T)
names(DESeq2)[1] <- "gene_name"

#Load and prepare shared OCR between two populations:
socr <- read.table(paste0("data/", list.files(path=paste0(getwd(),"/data"), pattern= "OCR")),
                   sep = "\t", dec = ".",header = TRUE, quote = "", stringsAsFactors = F)
socr <- socr[, "Region.ID", drop = F]
socr$Chr <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[1]))
socr$ranges <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[2]))
grocr <- GRanges(seqnames = socr$Chr, 
               ranges = socr$ranges, 
               strand = NULL)
#Create an empty dataframe to be filled with data over the loop
Overall_summary <- data.frame(matrix(ncol=3))
names(Overall_summary) <- c("Analysis",pop[2],pop[1])

for (i in c(1:length(pop))) {
  for (m in c("ChIP", "ATAC")){
  file_list <- list.files(path=paste0(getwd(),"/data"), pattern= m)
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
  Overall_summary[1,1] <- "ATAC_Overlapping_peaks_with_H3K27ac_ChIP"
  Overall_summary[1,4-i] <- nrow(olpeaks)
  
  #Obtain active enhancers by filtering overlapping peaks in promoters:
  prom <- promoters(TxDb)
  inpromoters <- findOverlaps(prom, gr3)
  act.enh <- olpeaks[-c(unique(subjectHits(inpromoters))), 1:3]
  Overall_summary[2,1] <- "Active_enhancers"
  Overall_summary[2,4-i] <- nrow(act.enh)
  
  #Export files in desired formats
  formats <- c(".txt", ".bed")
  col_names <- c(T,F)
  for (o in 1:length(formats)){
  write.table(act.enh, file = paste0("output/", pop[i] ,"_Active_enhancers",formats[o]),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  #Methylation in active enhancers
  gr4 <- GRanges(seqnames = act.enh$Chr, 
                 ranges = paste0(act.enh$Start,"-",act.enh$End), 
                 strand = NULL)
  gr5 <- GRanges(seqnames = DMR$Chr, 
                 ranges = paste0(DMR$Start,"-",DMR$End), 
                 strand = NULL,
                 Met.smp = DMR[,pop[2]],
                 Met.ctl = DMR[,pop[1]])
  overlap <- findOverlaps(gr5, gr4)
  act.enh.DMR_ <- act.enh[unique(subjectHits(overlap)),]
  write.table(act.enh.DMR_, file = paste0("output/", pop[i] ,"_Active_enhancers_not_promoter_with_DMR",formats[o]),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = col_names[o])
  Overall_summary[3,1] <- "Active_enhancers_with_DMR"
  Overall_summary[3,4-i] <- nrow(act.enh.DMR_)
  
  #do the overlap in the other direction
  overlap2 <- findOverlaps(gr4, gr5)
  DMR.act.enh <- DMR[unique(subjectHits(overlap2)),]
  write.table(DMR.act.enh, paste0("output/", pop[i] ,"_DMR_Overlapping_Active_enhancers",formats[o]),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = col_names[o])
  Overall_summary[6,1] <- "DMR_Overlapping_Active_enhancers"
  Overall_summary[6,4-i] <- nrow(DMR.act.enh)
  Overall_summary[7,1] <- "of_which_hypomethylated"
  if (i == 1){
  Overall_summary[7,4-i] <- nrow(DMR.act.enh[which(DMR.act.enh[,pop[2]]>DMR.act.enh[,pop[1]]),])
  } else {
    Overall_summary[7,4-i] <- nrow(DMR.act.enh[which(DMR.act.enh[,pop[2]]<DMR.act.enh[,pop[1]]),])
  }
  
  #Methylation in H3k27ac
  #Find overlapping peaks (population-specific H3K27ac mark + open region)
  overlap3 <- findOverlaps(grchip, grocr)
  openH3K27ac <- socr[unique(subjectHits(overlap3)),]
  grH3 <- GRanges(seqnames = openH3K27ac$Chr, 
                 ranges = openH3K27ac$ranges, 
                 strand = NULL)
  openH3K27ac <- separate(openH3K27ac, col = "ranges", into = c("Start", "End"), sep = "-", remove = T)
  openH3K27ac <- openH3K27ac[, c("Chr", "Start", "End")]
  write.table(openH3K27ac, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac",formats[o]),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = col_names[o])
  Overall_summary[8,1] <- "Shared_ATAC_H3K27ac"
  Overall_summary[8,4-i] <- nrow(openH3K27ac)
  
  #Obtain active enhancers by filtering overlapping peaks in promoters:
  inpromoters <- findOverlaps(prom, grH3)
  openH3K27acp <- openH3K27ac[-c(unique(subjectHits(inpromoters))), 1:3]
  write.table(openH3K27acp, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac_not_promoter",formats[o]),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = col_names[o])
  Overall_summary[9,1] <- "Shared_ATAC_H3K27ac_not_promoter"
  Overall_summary[9,4-i] <- nrow(openH3K27acp)
  
  overlap5 <- findOverlaps(gr5, grH3)
  H3K27open.DMR_ <- openH3K27ac[unique(subjectHits(overlap5)),]
  H3K27open.DMR_ <- H3K27open.DMR_[which(H3K27open.DMR_$Chr!="NA"),]
  write.table(H3K27open.DMR_, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac_not_promoter_with_DMR",formats[o]),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = col_names[o])
  Overall_summary[10,1] <- "Shared_ATAC_H3K27ac_not_promoter_with_DMR"
  Overall_summary[10,4-i] <- nrow(H3K27open.DMR_)
  
  overlap6 <- findOverlaps(grH3, gr5)
  DMR.H3K27open <- DMR[unique(subjectHits(overlap6)),]
  write.table(DMR.H3K27open, file = paste0("output/", pop[i] ,"_DMR_Overlapping_shared_ATAC_H3K27ac_not_promoter",formats[o]),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  Overall_summary[13,1] <- "DMR_Overlapping_shared_ATAC_H3K27ac_not_promoter"
  Overall_summary[13,4-i] <- nrow(DMR.H3K27open)
  Overall_summary[14,1] <- "s of_which_hypomethylated"
  if (i == 1){
    Overall_summary[14,4-i] <- nrow(DMR.act.enh[which(DMR.H3K27open[,pop[2]]>DMR.H3K27open[,pop[1]]),])
  } else {
    Overall_summary[14,4-i] <- nrow(DMR.act.enh[which(DMR.H3K27open[,pop[2]]<DMR.H3K27open[,pop[1]]),])
    }
  }
  
#Look for genes 100 kb around active enhancers with or without DMR using this loop
#Call files containing these regions from environtment as they have been previously created

aedmr <- grep("\\.DMR_", names(.GlobalEnv),value=TRUE)
mapped_genes <- mappedkeys(org.SYMBOL2EG)
xx <- as.list(org.SYMBOL2EG[mapped_genes])
bb <-as.data.frame(unlist(xx))
bb$gene_name <- row.names(bb)
names(bb) <- c("gene_id","gene_name")
for (ae in 1:length(aedmr)){
  prgr <- eval(as.symbol(grep("\\.DMR_", names(.GlobalEnv),value=TRUE)[ae]))
  gr8 <- GRanges(seqnames = prgr$Chr, 
                 ranges = paste0(prgr$Start,"-",prgr$End), 
                 strand = NULL)
  g.around.df <- data.frame(matrix(ncol = 6, nrow = 0))
  names(g.around.df) <- c("act.enh","Chr","Start", "End", "gene_strand", "gene_name")
  g.around.fdf <- g.around.df
  for (a in 1:length(gr8)){
    genes.around <- vector()
    gr100kb <- resize(gr8[a], width(gr)+100000, fix = "center")
    gene.by.act.enh <- subsetByOverlaps(genes(TxDb), gr100kb) 
    if (length(gene.by.act.enh)>0){
      df <- data.frame(act.enh = rep(a, length(gene.by.act.enh)),
                       Chr = rep(seqnames(gr8[a]), length(gene.by.act.enh)),
                       Start = rep(start(ranges(gr8[a])),length(gene.by.act.enh)),
                       End = rep(end(ranges(gr8[a])),length(gene.by.act.enh)),
                       gene_strand = strand(gene.by.act.enh),
                       gene_id = mcols(gene.by.act.enh))
      df$gene_id <- as.numeric(df$gene_id)
      df1 <- merge(x= df, y= bb, by="gene_id", all.x = T)
      df2 <- ddply(df1, .(act.enh), summarize, 
                   act.enh=paste(unique(act.enh),collapse=","),
                   Chr = paste(unique(Chr),collapse=","),
                   Start = paste(unique(Start),collapse=","),
                   End = paste(unique(End),collapse=","),
                   gene_strand = paste(gene_strand,collapse=","),
                   gene_name = paste(unique(gene_name),collapse=","))
      g.around.df <- rbind(g.around.df, df2)
      g.around.fdf <- rbind(g.around.fdf,df1)
    } else {
      df0 <-data.frame(act.enh = a,
            Chr = seqnames(gr8[a]),
            Start = start(ranges(gr8[a])),
            End = end(ranges(gr8[a])),
            gene_strand = "NA",
            gene_name = "No genes found")
      g.around.df <- rbind(g.around.df, df0)
     }
  }
  write.table(g.around.df, paste0("output/",pop[i], "_", aedmr[ae], "genes_around_100kb.txt"),
              sep = "\t", row.names = T, col.names = F, quote = F)
  #generate also excel file
  write_xlsx(g.around.df, paste0("output/",pop[i], "_", aedmr[ae], "genes_around_100kb.xlsx"))
 
  #Analysis of genes associated to active enhancers and their transcriptomic activity
  #careful! if names changes order may change
  trans_DMR <- merge(DESeq2, g.around.fdf, by.y = "gene_name", all.x = F)
  trans_DMR <- unique(trans_DMR[,c(1:8)])
  
  if (str_detect(aedmr[ae], "H3K27open.")){
    n <- 11 } else {
      n <- 4
    }
  if (i == 1){
    Overall_summary[n,3] <- nrow(subset(trans_DMR, trans_DMR$log2FoldChange>=2 & trans_DMR$padj<=0.01))
    Overall_summary[n+1,3] <- nrow(subset(trans_DMR, trans_DMR$log2FoldChange<=-2 & trans_DMR$padj<=0.01))
  } else {
    Overall_summary[n,2] <- nrow(subset(trans_DMR, trans_DMR$log2FoldChange>=2 & trans_DMR$padj<=0.01))
    Overall_summary[n+1,2] <- nrow(subset(trans_DMR, trans_DMR$log2FoldChange<=-2 & trans_DMR$padj<=0.01))
  }
  
  }
}

#Complete Overall_summary table
Overall_summary[4,1] <- paste("genes with log2FC>=2", pop[2], "vs", pop[1])
Overall_summary[5,1] <- paste("genes with log2FC<=-2", pop[2], "vs", pop[1])
Overall_summary[11,1] <- paste("s genes with log2FC>=2", pop[2], "vs", pop[1])
Overall_summary[12,1] <- paste("s genes with log2FC<=-2", pop[2], "vs", pop[1])

write_xlsx(Overall_summary, "output/Overall_summary_active_enhancers.xlsx")


#Do graph bar plots for summary
dfplot <- data.frame(matrix(ncol = 1, nrow= nrow(Overall_summary)*2))
dfplot$analysis <- Overall_summary$Analysis
dfplot$analysis <- as.character(dfplot$analysis)
dfplot[c(1:nrow(Overall_summary)),1] <- Overall_summary[,pop[2]]
dfplot[c(1:nrow(Overall_summary)),3] <- pop[2]
dfplot[c((nrow(Overall_summary)+1):(nrow(Overall_summary)*2)),1] <- Overall_summary[,pop[1]]
dfplot[c((nrow(Overall_summary)+1):(nrow(Overall_summary)*2)),3] <- pop[1]
names(dfplot) <- c("Hits","analysis","type")

hitsbar <- ggplot(dfplot, aes(y=analysis, x=Hits, fill=type)) + geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Hits), vjust=0.5, hjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5)+
  ylab("Type of analysis") + xlab("Num of hits(log10)")+ scale_x_log10(expand = expansion(add= c(0 , 1)))+
  # theme(axis.text.x = element_text(angle = 60, size = 10, hjust =1, face="bold"))+
  scale_y_discrete(limits = rev(Overall_summary$Analysis))+
  #scale_fill_manual(values = c("firebrick","gray"))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.3, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"))
#Export pdf with table and bar graph
pdf(file = "figs/Overall_summary.pdf", width = 10, height = 6)
grid.table(Overall_summary, rows = NULL)
print(hitsbar)
dev.off()
