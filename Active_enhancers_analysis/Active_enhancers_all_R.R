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

# Generate database for species to be studied
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
BM <- getBM (attributes=c("entrezgene_id", "chromosome_name", "start_position", "end_position" , "ensembl_gene_id", "external_gene_name", "strand"),
             mart = ensembl, verbose = T)
prom <- promoters(TxDb)

# Create Granges object from BM database
# Take out mithocondrial and weird chromosome annotations
BMgr <- subset(BM, as.numeric(chromosome_name)>= 1 & as.numeric(chromosome_name) <=50 | chromosome_name == "X" | chromosome_name == "Y")
# Nomenclature for chromosomes should be "chrX"
BMgr$chromosome_name <- sapply(BMgr$chromosome_name, function(x) {paste0("chr",x)})
BMgr$strand <- sapply(BMgr$strand, function(x){if (x ==1){print("+")} else {print("-")}})
# GRanges generated just in case, but not actually needed as we need bed file for bedtools window tool.
grBM <- GRanges(seqnames = BMgr$chromosome_name, 
                ranges = paste0 (BMgr$start_position, "-", BMgr$end_position), 
                strand = BMgr$strand,
                gene_name = BMgr$external_gene_name)
BMgenes <- BMgr[, c(2,3,4,6,7)]
names(BMgenes)[c(1:4)] <- c("Chr", "Start", "End", "gene_name")

# Load all files needed
# Load DMR file between two samples:
DMR <- read.table("data/DMR.txt",
                  sep = "\t", quote = "",
                  dec = ".", header = T, na.strings = T)
# Fix problem with chrchr9, usual issue?
DMR$Chr <- sapply(strsplit(as.character(DMR$Chr), split="chr", fixed=TRUE), function(x){print(x[2])})
DMR$Chr <- sapply(DMR$Chr, function(x){if (x == "") {print("chr9")} else {paste0("chr", x)}})
names(DMR) <- c("Chr", "Start", "End", pop[1], pop[2])
DMR <- DMR[order(as.numeric(gsub("chr", "", DMR$Chr)), 
                 as.numeric(DMR$Start),
                 decreasing = F, na.last = T), ]

# Add to your data directory the DESeq2 file comparing your samples to analyse
DESeq2 <- read.table (file = paste0("data/", list.files(path=paste0(getwd(),"/data"),
                                                        pattern= "DESeq2", ignore.case = T)),
                      sep = "\t", quote = "", dec = ".", header=T)
names(DESeq2)[1] <- "gene_name"

# Load and prepare shared OCR between two populations:
socr <- read.table(paste0("data/", list.files(path=paste0(getwd(),"/data"), pattern= "Shared", ignore.case = T)),
                   sep = "\t", dec = ".", header = TRUE, quote = "", stringsAsFactors = F)
socr <- socr[, "Region.ID", drop = F]
socr$Chr <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[1]))
socr$ranges <- sapply(strsplit(socr$Region.ID, split=':', fixed=TRUE), function(x) (x[2]))
grocr <- GRanges(seqnames = socr$Chr, 
                 ranges = socr$ranges, 
                 strand = NULL)

# Create an empty dataframe to be filled with data over the loop
Overall_summary <- data.frame(matrix(ncol=3))
names(Overall_summary) <- c("Analysis",pop[2],pop[1])

# Create function to look for genes around elements such as active enhancers
look.around <- function(x){
  g_around_more_rows <- data.frame(matrix(ncol = 6, nrow = 0))
  names(g_around_more_rows) <- c("element","Chr","Start", "End", "gene_strand", "gene_name")
  g_around <- g_around_more_rows
  for (a in 1:length(x)){
    gr100kb <- resize(x[a], width(x[a])+100000, fix = "center")
    gene.by.act.enh <- subsetByOverlaps(grBM, gr100kb)
    if (length(gene.by.act.enh)>0){
      df1 <- data.frame(element = rep(a, length(gene.by.act.enh)),
                        Chr = rep(seqnames(x[a]), length(gene.by.act.enh)),
                        Start = rep(start(ranges(x[a])), length(gene.by.act.enh)),
                        End = rep(end(ranges(x[a])), length(gene.by.act.enh)),
                        gene_strand = strand(gene.by.act.enh),
                        gene_name = mcols(gene.by.act.enh))
      df2 <- ddply(df1, .(element), summarize, 
                   element=paste(unique(element), collapse=","),
                   Chr = paste(unique(Chr), collapse=","),
                   Start = paste(unique(Start), collapse=","),
                   End = paste(unique(End), collapse=","),
                   gene_strand = paste(gene_strand, collapse=","),
                   gene_name = paste(unique(gene_name), collapse=","))
    } else {
      df1 <-data.frame(element = a,
                       Chr = seqnames(x[a]),
                       Start = start(ranges(x[a])),
                       End = end(ranges(x[a])),
                       gene_strand = "NA",
                       gene_name = "No genes found")
      df2 <- df1
    }
    g_around_more_rows <- rbind(g_around_more_rows,df1)
    g_around <- rbind(g_around, df2)
  }
  return(list(g_around_more_rows,g_around))
}

for (i in c(1:length(pop))) {
  for (m in c("ChIP", "ATAC")){
    file_list <- list.files(path=paste0(getwd(),"/data"), pattern= m)
    # Read file tables
    tble <- read.table(paste0("data/", grep(pop[i], file_list, value = T)), sep = "\t", quote = "",
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
                    ranges = paste0(tble$Start, "-", tble$End), 
                    strand = NULL,
                    `-log10.pval`= tble$`-log10.pval`,
                    FoldEnrichment = tble$FoldEnrichment,
                    Anno.Gene = tble$Anno.Gene,
                    Peak.location = tble$Gene.section,
                    Distance.to.TSS = tble$Distance.to.TSS)
      grchip <- GRanges(seqnames = tble$Chr, 
                        ranges = paste0(tble$Start, "-", tble$End), 
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
                    ranges = paste0(atac$Start, "-", atac$End), 
                    strand = NULL,
                    Anno.Gene = atac$Anno.Gene)
    }
    assign(paste0(m,".", pop[i], ".gr"), gr)
  }
  
  # Find overlapping peaks
  # CAREFUL: order of objects in `findOverlaps` matters!
  overlap <- findOverlaps(eval(as.symbol(grep(paste0("ChIP.", pop[i]), names(.GlobalEnv), value=TRUE))), 
                          eval(as.symbol(grep(paste0("ATAC.", pop[i]), names(.GlobalEnv), value=TRUE))))
  olpeaks <- atac[unique(subjectHits(overlap)),]
  gr3 <- GRanges(seqnames = olpeaks$Chr, 
                 ranges = paste0(olpeaks$Start, "-", olpeaks$End), 
                 strand = NULL)
  Overall_summary[1,1] <- "ATAC_Overlapping_peaks_with_H3K27ac_ChIP"
  Overall_summary[1,4-i] <- nrow(olpeaks)
  
  # Obtain active enhancers by filtering overlapping peaks in promoters:
  inpromoters <- findOverlaps(prom, gr3)
  act.enh <- olpeaks[-c(unique(subjectHits(inpromoters))), 1:3]
  Overall_summary[2,1] <- "Active_enhancers"
  Overall_summary[2,4-i] <- nrow(act.enh)
  
  # Annotate already in GRanges objects with active enhancers (without promoters):
  gr4 <- GRanges(seqnames = act.enh$Chr, 
                 ranges = paste0(act.enh$Start, "-", act.enh$End), 
                 strand = NULL)
  g_r <- look.around(gr4)
  g_around_more_rows <- as.data.frame(g_r[1]) 
  g_around <- as.data.frame(g_r[2])
  
  # Export files in .txt format as bed will not be needed in this script
  write.table(g_around, 
              file = paste0("output/", pop[i] ,"_Active_enhancers_annotated", ".txt"),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  
  # Methylation in active enhancers
  gr5 <- GRanges(seqnames = DMR$Chr, 
                 ranges = paste0(DMR$Start, "-", DMR$End), 
                 strand = NULL,
                 Met.smp = DMR[,pop[2]],
                 Met.ctl = DMR[,pop[1]])
  overlap <- findOverlaps(gr5, gr4)
  act.enh.DMR <- g_around[unique(subjectHits(overlap)),]
  act.enh.DMR.mr <- merge(act.enh.DMR, g_around_more_rows, all.y = T)
  write.table(act.enh.DMR, file = paste0("output/", pop[i] ,"_Active_enhancers_with_DMR_annotated", ".txt"),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  Overall_summary[3,1] <- "Active_enhancers_with_DMR"
  Overall_summary[3,4-i] <- nrow(act.enh.DMR)
  
  # Do the overlap in the other direction
  overlap2 <- findOverlaps(gr4, gr5)
  DMR.act.enh <- DMR[unique(subjectHits(overlap2)),]
  grdmr <- GRanges(seqnames = DMR.act.enh$Chr, 
                   ranges = paste0(DMR.act.enh$Start, "-", DMR.act.enh$End), 
                   strand = NULL)
  g_r <- look.around(grdmr)
  dmr_g_around_more_rows <- as.data.frame(g_r[1])
  dmr_g_around <- as.data.frame(g_r[2])
  dmr_g_aroundc <- merge(dmr_g_around, DMR.act.enh)
  
  write.table(dmr_g_aroundc, paste0("output/", pop[i] ,"_DMR_Overlapping_Active_enhancers_annotated",".txt"),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  Overall_summary[6,1] <- "DMR_Overlapping_Active_enhancers"
  Overall_summary[6,4-i] <- nrow(DMR.act.enh)
  Overall_summary[7,1] <- "of_which_hypomethylated"
  if (i == 1){
    Overall_summary[7,4-i] <- nrow(DMR.act.enh[which(DMR.act.enh[,pop[2]]>DMR.act.enh[,pop[1]]),])
  } else {
    Overall_summary[7,4-i] <- nrow(DMR.act.enh[which(DMR.act.enh[,pop[2]]<DMR.act.enh[,pop[1]]),])
  }
  
  # Methylation in H3k27ac
  # Find overlapping peaks (population-specific H3K27ac mark + open region)
  overlap3 <- findOverlaps(grchip, grocr)
  openH3K27ac <- socr[unique(subjectHits(overlap3)),]
  grH3 <- GRanges(seqnames = openH3K27ac$Chr, 
                  ranges = openH3K27ac$ranges, 
                  strand = NULL)
  Overall_summary[8,1] <- "Shared_ATAC_H3K27ac"
  Overall_summary[8,4-i] <- nrow(openH3K27ac)
  
  # Obtain active enhancers by filtering overlapping peaks in promoters:
  inpromoters <- findOverlaps(prom, grH3)
  openH3K27acp <- openH3K27ac[-c(unique(subjectHits(inpromoters))), 1:3]
  grH3wop <- GRanges(seqnames = openH3K27acp$Chr, 
                     ranges = openH3K27acp$ranges, 
                     strand = NULL)
  g_r <- look.around(grH3wop)
  H3_g_around_more_rows <- as.data.frame(g_r[1])
  H3_g_around <- as.data.frame(g_r[2])
  
  write.table(H3_g_around, file = paste0("output/", pop[i] ,"_shared_ATAC_H3K27ac_not_promoter_annotated", ".txt"),
              sep = "\t", quote = F, dec = ".", row.names = F, col.names = T)
  Overall_summary[9,1] <- "Shared_ATAC_H3K27ac_not_promoter"
  Overall_summary[9,4-i] <- nrow(openH3K27acp)
  
  overlap5 <- findOverlaps(gr5, grH3wop)
  H3K27open.DMR <- H3_g_around[unique(subjectHits(overlap5)),]
  H3K27open.DMR.mr <- merge(H3K27open.DMR, H3_g_around_more_rows, all.y = T)
  write.table(H3K27open.DMR, file = paste0("output/annotation/", pop[i],
                                           "_shared_ATAC_H3K27ac_not_promoter_with_DMR_annotated", ".txt"),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  Overall_summary[10,1] <- "Shared_ATAC_H3K27ac_not_promoter_with_DMR"
  Overall_summary[10,4-i] <- nrow(H3K27open.DMR)
  
  overlap6 <- findOverlaps(grH3wop, gr5)
  DMR.H3K27open <- DMR[unique(subjectHits(overlap6)),]
  grdmrH3 <- GRanges(seqnames = DMR.H3K27open$Chr, 
                     ranges = paste0(DMR.H3K27open$Start, "-", DMR.H3K27open$End), 
                     strand = NULL)
  g_r <- look.around(grdmrH3)
  dmrH3_g_around_more_rows <- as.data.frame(g_r[1])
  dmrH3_g_around <- as.data.frame(g_r[2])
  dmrH3_g_aroundc <- merge(dmrH3_g_around, DMR.H3K27open) # if willing to include all rows, just merge with *_more_rows (all.y = T)
  
  write.table(dmrH3_g_aroundc, file = paste0("output/", pop[i] ,"_DMR_Overlapping_shared_ATAC_H3K27ac_not_promoter_annotated",".txt"),
              sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)
  
  Overall_summary[13,1] <- "DMR_Overlapping_shared_ATAC_H3K27ac_not_promoter"
  Overall_summary[13,4-i] <- nrow(DMR.H3K27open)
  Overall_summary[14,1] <- "s of_which_hypomethylated"
  if (i == 1){
    Overall_summary[14,4-i] <- nrow(DMR.act.enh[which(DMR.H3K27open[,pop[2]]>DMR.H3K27open[,pop[1]]),])
  } else {
    Overall_summary[14,4-i] <- nrow(DMR.act.enh[which(DMR.H3K27open[,pop[2]]<DMR.H3K27open[,pop[1]]),])
  }
  
  # Analysis of genes associated to active enhancers and their transcriptomic activity
  # CAREFUL! If names change, order may change
  
  for (mr in 1:length(grep("\\.mr", names(.GlobalEnv), value=TRUE))){
    trans_DMR <- merge(DESeq2, eval(as.symbol(grep("\\.mr", names(.GlobalEnv), value=TRUE)[mr])), by.y = "gene_name", all.x = F)
    trans_DMR <- unique(trans_DMR[,c(1:8)])
    
    if (str_detect(grep("\\.mr", names(.GlobalEnv), value=TRUE)[mr], "H3K27open.")){
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

# Complete Overall_summary table
Overall_summary[4,1] <- paste("genes with log2FC>=2", pop[2], "vs", pop[1])
Overall_summary[5,1] <- paste("genes with log2FC<=-2", pop[2], "vs", pop[1])
Overall_summary[11,1] <- paste("s genes with log2FC>=2", pop[2], "vs", pop[1])
Overall_summary[12,1] <- paste("s genes with log2FC<=-2", pop[2], "vs", pop[1])

write_xlsx(Overall_summary, "output/Overall_summary_active_enhancers_all_R.xlsx")


# Make graph (bar plot) for summary
dfplot <- data.frame(matrix(ncol = 1, nrow= nrow(Overall_summary)*2))
dfplot$analysis <- Overall_summary$Analysis
dfplot$analysis <- as.character(dfplot$analysis)
dfplot[c(1:nrow(Overall_summary)),1] <- Overall_summary[,pop[2]]
dfplot[c(1:nrow(Overall_summary)),3] <- pop[2]
dfplot[c((nrow(Overall_summary)+1):(nrow(Overall_summary)*2)),1] <- Overall_summary[,pop[1]]
dfplot[c((nrow(Overall_summary)+1):(nrow(Overall_summary)*2)),3] <- pop[1]
names(dfplot) <- c("Hits","analysis","type")

hitsbar <- ggplot(dfplot, aes(y=analysis, x=Hits, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Hits), vjust=0.5, hjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5) +
  ylab("Type of analysis") + xlab("Num of hits(log10)") + 
  scale_x_log10(expand = expansion(add= c(0 , 1))) +
  # theme(axis.text.x = element_text(angle = 60, size = 10, hjust =1, face="bold")) +
  scale_y_discrete(limits = rev(Overall_summary$Analysis)) +
  # scale_fill_manual(values = c("firebrick","gray")) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "grey",
                                        size = 0.3, linetype = "solid"),
        panel.grid.major = element_line(size = 0.05, linetype = 'solid',
                                        colour = "grey"))
# Export pdf with table and bar graph
pdf(file = "figs/Overall_summary.pdf", width = 10, height = 6)
grid.table(Overall_summary, rows = NULL)
print(hitsbar)
dev.off()
