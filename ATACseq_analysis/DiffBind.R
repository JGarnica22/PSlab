#This script is to perform DiffBind for differential analysis of peaks
#It was created for R version 4.0.2 2020-10-15

# Check if required packages are installed, if not install:
cran.packages <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "devtools", "rlang", "dplyr")
for (i in cran.packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    print(paste(i,"just installed"))
  }
  else {
    print(paste(i,"was already installed"))
  }
}

bioc.packages <- c("DESeq2", "biomaRt", "DiffBind","GenomicRanges", "Gviz", "trackViewer", "rtracklayer", 
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

# Load packages:
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(biomaRt)
library(DiffBind)
library(GenomicRanges)
library(trackViewer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(plyr)
library(dplyr)

# Set your working directory (the project you are working in):
setwd("C:/Users/jgarn/OneDrive - Universitat de Barcelona/Documentos/Bioinformatics/ATACseq/Diffbind")

## Indicate populations of interest, to be compared:
# First indicate control population, then sample:
pop <- c("Tconv", "Tet")

# Indicate species
species <- "mouse"

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

# Generate peakset data, you need bed files as well as bam and respective bai files.
list_files <- list.files(path= paste0(getwd(), "/data/"),
                       pattern= "*.bed")
list_files <- sapply(c(1:length(list_files)), function(x){strsplit(as.character(list_files), "_")[[x]][2]})

samples <- data.frame(matrix(nrow=length(list_files), ncol=8))
names(samples) <- c("SampleID", "Condition", "Factor", "Treatment", "Replicate", "bamReads", "Peaks",
                    "PeakCaller")
samples$SampleID <- list_files
samples$Condition <- gsub('.{1}$', '', list_files)
samples$Factor <- samples$Condition
samples$Treatment <- samples$Condition
samples[c(grep(pop[1], list_files)),"Replicate"] <- 1:length(grep(pop[1], list_files))
samples[c(grep(pop[2], list_files)),"Replicate"] <- 1:length(grep(pop[2], list_files))
samples$PeakCaller <- "bed"
list_bed <- list.files(path= paste0(getwd(), "/data"),
                         pattern= "*.bed")
samples$Peaks <-sapply("data/", paste0, list_bed)
list_bam <- list.files(path= paste0(getwd(), "/data"),
                       pattern= "*.bam$")
samples$bamReads <-sapply("data/", paste0, list_bam)

dbdata <- dba(sampleSheet=samples)
plot(dbdata)
dba.plotPCA(dbdata, label = "ID")

# Counting reads
# Next step is to calculate a binding matrix based on the read counts for every sample rather than only based on the
# peaks called. Let's generate the counts and have a look at the data:

dbdata.count <- dba.count(dbdata, summits = 250)
plot(dbdata.count)
dba.plotPCA(dbdata.count, label = "ID")


# Establishing a contrast
# Before running the differential analysis, we need to tell DiffBind which cell lines fall in which groups. 
dbdata.contrast <- dba.contrast(dbdata.count, categories=DBA_CONDITION, minMembers = 2)

# Differential analysis
dbdata.anal <- dba.analyze(dbdata.contrast, method=DBA_ALL_METHODS)
dba.plotPCA(dbdata.anal, contrast=1, method=DBA_ALL_METHODS, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbdata.anal, contrast=1, method=DBA_ALL_METHODS)
dba.plotHeatmap(dbdata.anal, method=DBA_ALL_METHODS)


# Retrieve differentially bound sites

dbdata.DB <- dba.report(dbdata.anal, method=DBA_ALL_METHODS)

# Look for genes around these regions

BM <- getBM(attributes=c("entrezgene_id", "chromosome_name", "start_position", "end_position" , "ensembl_gene_id", "external_gene_name", "strand"),
             mart = ensembl, verbose = T)
BMgr <- subset(BM, as.numeric(chromosome_name)>= 1 & as.numeric(chromosome_name) <=50 | chromosome_name == "X" | chromosome_name == "Y")
# Nomenclature for chromosomes should be "chrX"
#BMgr$chromosome_name <- sapply(BMgr$chromosome_name, function(x) {paste0("chr",x)})
BMgr$strand <- sapply(BMgr$strand, function(x){if (x ==1){print("+")} else {print("-")}})
# GRanges generated just in case, but not actually needed as we need bed file for bedtools window tool.
grBM <- GRanges(seqnames = BMgr$chromosome_name, 
                ranges = paste0 (BMgr$start_position, "-", BMgr$end_position), 
                strand = BMgr$strand,
                gene_name = BMgr$external_gene_name)

look.around <- function(x){
  g_around_more_rows <- data.frame(matrix(ncol = 6, nrow = 0))
  names(g_around_more_rows) <- c("element","Chr","Start", "End", "gene_strand", "gene_name")
  g_around <- g_around_more_rows
  for (a in 1:length(x)){
    gr100kb <- resize(x[a], width(x[a])+5000, fix = "center")
    gene.by.act.enh <- subsetByOverlaps(grBM, gr100kb)
    if (length(gene.by.act.enh)>0){
      df1 <- data.frame(element = rep(a, length(gene.by.act.enh)),
                        Chr = rep(seqnames(x[a]), length(gene.by.act.enh)),
                        Start = rep(start(ranges(x[a])), length(gene.by.act.enh)),
                        End = rep(end(ranges(x[a])), length(gene.by.act.enh)),
                        gene_strand = strand(gene.by.act.enh),
                        gene_name = mcols(gene.by.act.enh),
                        Fold=mcols(x[a])[4],
                        p.value=mcols(x[a])[5],
                        FDR=mcols(x[a])[6])
      df2 <- ddply(df1, .(element), summarize, 
                   element=paste(unique(element), collapse=","),
                   Chr = paste(unique(Chr), collapse=","),
                   Start = paste(unique(Start), collapse=","),
                   End = paste(unique(End), collapse=","),
                   gene_strand = paste(gene_strand, collapse=","),
                   gene_name = paste(unique(gene_name), collapse=","),
                   Fold=paste(unique(Fold), collapse=","),
                   p.value=paste(unique(p.value), collapse=","),
                   FDR=paste(unique(FDR), collapse=",")
                   )
    } else {
      df1 <-data.frame(element = a,
                       Chr = seqnames(x[a]),
                       Start = start(ranges(x[a])),
                       End = end(ranges(x[a])),
                       gene_strand = "NA",
                       gene_name = "No genes found",
                       Fold=mcols(x[a])[4],
                       p.value=mcols(x[a])[5],
                       FDR=mcols(x[a])[6])
      df2 <- df1
    }
    g_around_more_rows <- rbind(g_around_more_rows,df1)
    g_around <- rbind(g_around, df2)
    }
  return(list(g_around_more_rows,g_around))
}

look.around(dbdata.DB)

