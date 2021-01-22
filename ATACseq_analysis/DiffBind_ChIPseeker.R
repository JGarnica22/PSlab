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
                   "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "ggplot", "ChIPseeker", "clusterProfiler",
                   "ChIPpeakAnno", "UpsetR", "ggupset", "ggimage", "ReactomePA")
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
library(tidyverse)
library(gridExtra)
library(edgeR)
library(ggthemes)
library(readxl)
library(writexl)
library(ChIPseeker)
library(clusterProfiler)
library(ChIPpeakAnno)
library(UpSetR)
library(ggupset)
library(ggimage)
library(ReactomePA)


# Set your working directory (the project you are working in):
setwd("~/Downloads/ATACseq_pipeline")

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

# Reading in the peaksets
# The easiest way to set up an experiment to analyze is with a sample sheet. The sample sheet
# can be a dataframe, or it can be read directly from a csv file.

# Generate automatically sample sheet peakset data, you need peak files as well as bam and respective bai files.
# Examine your files and generate short names for your samples based on the names of bam files
list_files <- list.files(path= paste0(getwd(), "/data"),
                       pattern= "*.bam$")
list_files <- sapply(c(1:length(list_files)), function(x){strsplit(as.character(list_files), "_")[[x]][2]})
# Aquí ha de ser el string [2] perquè agafi el nom de la població: "nd_Tconv1_24121_..."

samples <- data.frame(matrix(nrow=length(list_files), ncol=8))
names(samples) <- c("SampleID", "Condition", "Factor", "Treatment", "Replicate", "bamReads", "Peaks",
                    "PeakCaller")
samples$SampleID <- list_files
samples$Condition <- gsub('.{1}$', '', list_files)
# You may change factor and treatment if experiment presented different subgroups.
samples$Factor <- samples$Condition
samples$Treatment <- samples$Condition
samples[c(grep(pop[1], list_files)),"Replicate"] <- 1:length(grep(pop[1], list_files))
samples[c(grep(pop[2], list_files)),"Replicate"] <- 1:length(grep(pop[2], list_files))
# Peakcaller will be changing depending on the file use (see DiffBind vignette)
# To call macs2 output format: 
# "bed": .bed file; peak score is in fifth column; 
# "narrow": default peak.format: narrowPeaks file;
# "macs": MACS .xls file.
samples$PeakCaller <- "macs"
list_bed <- list.files(path= paste0(getwd(), "/data"),
                         pattern= "*peaks.xls")
samples$Peaks <-as.vector(sapply("data/", paste0, list_bed))
list_bam <- list.files(path= paste0(getwd(), "/data"),
                       pattern= "*.bam$")
samples$bamReads <-as.vector(sapply("data/", paste0, list_bam))


# Constructs a new DBA object from the sample sheet
# This contains how many peaks are in each peakset, as well as (in the first line) the total number
# of unique peaks after merging overlapping ones
dbdata <- dba(sampleSheet=samples)
# Using the data from the peak calls, a correlation heatmap and PCA can be generated which gives an
# initial clustering of the samples using the cross-correlations of each row of the binding matrix
pdf(file = paste0("figs/DiffBind_plots_peak_caller_score.pdf"))
plot(dbdata, main="Correlation heatmap peak caller score")
pca1 <- dba.plotPCA(dbdata, label = "ID")
# Also we can compare the number of peaks (intervals) in a bar plot
dba.show(dbdata) %>% ggplot(aes(x=ID, y=Intervals)) + geom_bar(stat="identity") +
ggtitle("Number of peaks") + xlab("Sample") + theme_economist()
# and a Venn graph showing the binding site overlaps between all samples:
dba.plotVenn(dbdata, dbdata$masks$All)
grid.newpage()
grid.table(dba.show(dbdata))
dev.off()

# Counting reads
# Next step is to calculate a binding matrix based on the read counts for every sample (affinity scores) rather than only based on the
# peaks called. Let's generate the counts and have a look at the data:
dbdata.count <- dba.count(dbdata, bUseSummarizeOverlaps=TRUE)
# This shows that all the samples are using the same number of consensus peakset. Also, a
# new column has been added, called FRiP, which stands for Fraction of Reads in Peaks. This
# is the proportion of reads for that sample that overlap a peak in the consensus peakset, and
# can be used to indicate which samples show more enrichment overall.

# We can also plot a new correlation graphs based on the affinity scores
pdf(file = paste0("figs/DiffBind_plots_affinity.pdf"))
plot(dbdata.count, main="Correlation heatmap with affinity data")
pca2 <- dba.plotPCA(dbdata.count, label = "ID")
grid.newpage()
grid.table(dba.show(dbdata.count))
dev.off()

# Normalizing the data:
norm <- dba.normalize(dbdata.count, normalize = DBA_NORM_LIB)

# Establishing a contrast
# Before running the differential analysis, we need to tell DiffBind which cell lines fall in which groups. 
dbdata.contrast <- dba.contrast(dbdata.count, categories=DBA_CONDITION, minMembers = 2 )
#Set threshold to use in analyse, let's get all peaks and then filter by pvalue or FDR
dbdata.contrast$config$th <- 1

# Differential analysis
# This will run an DESeq2 analysis using the default binding matrix
dbdata.anal <- dba.analyze(dbdata.contrast, method=DBA_ALL_METHODS)

# This shows how many sites are identified as being significantly differentially bound (DB) using the
# default threshold of FDR <= 0.05

# Unsing only differentially bound sites, we can see how our sample cluster
# IMPORTANT: This plot is not a "result" in the sense that the analysis is selecting for sites that differ between the two
# conditions, and hence are expected to form clusters representing the conditions.
pdf(file = paste0("figs/DiffBind_plots_differ_bound.pdf"))
dba.plotPCA(dbdata.anal, contrast=1, method=DBA_ALL_METHODS, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbdata.anal, contrast=1, method=DBA_ALL_METHODS, main="Significantly differentially bound sites")
grid.newpage()
grid.table(dba.show(dbdata.anal, bContrasts = T))
# We can also visualize how many differential bound sites are found based on the analysis performed
dba.plotVenn(dbdata.anal, contrast=1, method=DBA_ALL_METHODS)


# MA plots are a useful way to visualize the effect of normalization on data, as well as seeing
# which of the datapoints are being identified as differentially bound.
#  Each point represents a binding site, with points in red representing sites identified as differentially bound.
dba.plotMA(dbdata.anal, method=DBA_EDGER)
dba.plotMA(dbdata.anal, bXY=T, method = DBA_EDGER)
# Aquest segon plot no funciona

# Volcano plots also highlight significantly differentially bound sites and show their fold changes:
dba.plotVolcano(dbdata.anal, contrast=1, method = DBA_EDGER)
# Boxplots provide a way to view how read distributions differ between classes of binding sites
dba.plotBox(dbdata.anal)

# Binding affinity heatmap showing affinities for differentially bound sites
hmap <- colorRampPalette(c("red", "black", "green"))(n = 4)
dba.plotHeatmap(dbdata.count, method=DBA_ALL_METHODS, correlations=FALSE,
                                scale="row", colScheme = hmap, main="Affinities for differentially bound sites")
dev.off()

# Retrieve differentially bound sites
# These are returned as a GRanges object, appropriate for downstream processing.
# Here you can filter your data based on FDR threshold (th) and fold change (fold) or you can filter this later as dataframe

dbdata.DB <- dba.report(dbdata.anal) #, bUsePval = T, th = 0.05, fold = 2)

#write.table(as.data.frame(dbdata.DB), "out/DESeq2_peaks.txt", quote = F, sep = "\t", row.names = F)

## Represent results with ChIPseeker
# Prepare the TSS regions for your genome
# Determine upstream and downstream values
lim=5000
promoter <- getPromoters(TxDb=TxDb, upstream=lim, downstream=lim)
pdf("figs/Figures_peaks_seeker.pdf")
for (p in 1:length(list_bed)){
# Show where peaks fall on the chromosomes
e <- read.delim(paste0("data/",list_bed[p]), comment.char = "#") %>%  mutate(chr = sapply("chr", paste0, chr)) %>% 
toGRanges()
print(covplot(e, weightCol="pileup", title = paste0(strsplit(as.character(list_bed), "_")[[p]][5],
                                                           " peaks on chromosomes")))
## Profile of ChIP peaks binding to TSS regions
# Show how peaks fall around TSS regions
tagMatrix <- getTagMatrix(e, windows=promoter) 
tagHeatmap(tagMatrix, xlim=c(-lim, lim), color="red", 
           title = paste0(strsplit(as.character(list_bed), "_")[[p]][5]," peaks around TSS"), 
           xlab = "Genomic Region (5'->3')" )

# Average profile peaks binding to TSS region
plotAvgProf(tagMatrix, xlim=c(-lim, lim), conf = 0.95,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
  ggtitle(paste0(strsplit(as.character(list_bed), "_")[[p]][5]," average profile peaks around TSS"))
# Plots comparing differents samples
x <- read.delim(paste0("data/",list_bed[p]), comment.char = "#") %>%  mutate(chr = sapply("chr", paste0, chr)) %>% toGRanges()
assign(paste0("Granges_", strsplit(as.character(list_bed), "_")[[p]][5]), x)
}

# Error després del loop

po <- grep("Granges_*", names(.GlobalEnv),value=TRUE)
peaks <- GRangesList(Tconv1=eval(as.symbol(po[1])), Tconv3=eval(as.symbol(po[2])), Tet1=eval(as.symbol(po[3])), 
                     Tet3=eval(as.symbol(po[4])))
#Average profiles
col <- c("darkolivegreen", "darkolivegreen4", 
         "firebrick1", "darkred")
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-lim, lim)) + scale_color_manual(values= col) +
  ggtitle("All samples average profile peaks around TSS")
plotAvgProf(tagMatrixList, xlim=c(-lim, lim), conf=0.95,resample=500, facet="row") + scale_color_manual(values= col) +
  ggtitle("All samples average profile peaks around TSS")
tagHeatmap(tagMatrixList, xlim=c(-lim, lim), color=NULL, 
              title = names(tagMatrixList))

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=TxDb,
                       tssRegion=c(-lim, lim), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

# Compare functional profiles
# genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) <- sub("_", "\n", names(genes))
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            pvalueCutoff  = 0.05,
#                            pAdjustMethod = "BH")
# dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()

## Peak annotation with ChIPseeker, in this case let's annotate differentially bound peaks
# to add "chr" to seqnames use:
seqlevelsStyle(dbdata.DB) <- "UCSC"
peakAnno <- annotatePeak(dbdata.DB, tssRegion=c(-lim, lim),
                         TxDb=TxDb, annoDb="org.Mm.eg.db")
consensus <- as.data.frame(peakAnno)

#filter your data
th <- 0.05
FC <- 2
shared <- consensus %>% filter(p.value > th & abs(Fold) < FC)
Differ_Tconv <- consensus %>% filter(p.value <= th & Fold >= FC)
Differ_Tet <- consensus %>% filter(p.value <= th & Fold <= -FC)

#Export annotated data
write_xlsx(consensus, "out/All_consensus_peaks_annotated.xlsx")
write_xlsx(shared, "out/Shared_peaks_annotated.xlsx")
write_xlsx(Differ_Tconv, "out/Diff_Tconv_peaks_annotated.xlsx")
write_xlsx(Differ_Tet, "out/Diff_Tet_peaks_annotated.xlsx")


pathway <- enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
pdf("figs/Differentially_bound_annotation.pdf", width = 12)
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno, vennpie=T)
grid.newpage()
grid.table(as.data.frame(pathway))
dev.off()


# Homer annotation
# save report as bed file, regardless of the genome used, chrosome column must include "chr"
report <- as.data.frame(dbdata.DB)
# Export bed file for Homer annotatePeaks.pl to be run in linux:
report_bed <- report %>% mutate(Unique_ID=row.names(report)) %>% select(c(1:3, 12)) %>% 
  mutate(seqnames = sapply("chr", paste0, seqnames) )
report_bed[,5:6] <- NA
write.table (report_bed, "out/diffbind_report.bed",
             sep = "\t", dec = ".", quote = F, row.names = F, col.names = F)

## Run annotatePeaks from Homer on linux, this option also run gene ontology (GO) analysis
# annotatePeaks.pl diffbind_report.bed mm10 -go annotatepeak/GO > annotatepeak/homer_anno2.txt

#Import annotation results and merge with report dataframe
anno <- read.table("data/homer_anno2.txt",
                   header = T,
                   sep = "\t",
                   quote = "", 
                   dec = ".")
names(anno)[1] <- "PeakID"  
comp <- merge(mutate(report, PeakID=row.names(report)), anno, by = "PeakID", all = T)
comp <- comp[, -c(2:4, 6)]

#Export results
write_xlsx(comp, paste0("out/Diff_bound_annotated_", pop[1], "_v_", pop[2] ,".xlsx"))
