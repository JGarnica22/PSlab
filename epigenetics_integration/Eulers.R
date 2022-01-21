# Example for generation of eulers

#Libs

library(readxl)
library(xlsx)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(UpSetR)
library(ggplot2)
library(ggpubr)
library(gplots)
library(VennDiagram)
library(DiffBind)
library(ChIPpeakAnno)
library(seqsetvis)
library(gtools)
library(eulerr)


#Filter the differential peaks respect to TH0
tet_k4me3 <- t_k4me3[t_k4me3$FDR.step.up..Tet.vs..Th0. < 0.01 & t_k4me3$Fold.change..Tet.vs..Th0. > 0,] %>% 
  select(Region.ID) %>% separate(Region.ID, c("chr", "start", "end"), sep = ":|-")
tfh_k4me3 <- t_k4me3[t_k4me3$FDR.step.up..TFH.vs..Th0. < 0.01 & t_k4me3$Fold.change..TFH.vs..Th0. > 0,] %>% 
  select(Region.ID) %>% separate(Region.ID, c("chr", "start", "end"), sep = ":|-")

#Generate granges objects and perform the overlap
gr_tet <- GRanges(seqnames = tet_k4me3$chr,
            ranges = paste0(tet_k4me3$start, "-", tet_k4me3$end),
            strand = NULL)
gr_tfh <- GRanges(seqnames = tfh_k4me3$chr,
            ranges = paste0(tfh_k4me3$start, "-", tfh_k4me3$end),
            strand = NULL)

grs <- list(gr_tet, gr_tfh)
olaps <- ssvOverlapIntervalSets(grs)
colnames(mcols(olaps)) <- c("TET", "TFH")
ut <- mcols(olaps)
ut <- as.data.frame(ut)
for (k in 1:length(colnames(ut))){
  ut[,k] <- as.integer(as.logical(ut[,k]))
}

f <- euler(ut)

pdf("figs/olaps/k4me3_euleans_atac_19_01_2022.pdf", width = 15, height = 11)
plot(f,
     quantities = TRUE,
     fills = list(fill = c("#FF4848", "#FFD371", "#C2FFD9", "#9DDAC6")),
     labels = c("TET","TFH"),list(font = 4),
     main = list("Euler Diagram Overlaps K4me3 V1"))

plot(f,
     quantities = TRUE,
     fills = list(fill = c("#FCFFA6", "#C1FFD7", "#B5DEFF", "#CAB8FF")),
     lty = 1:4,
     labels = c("TET","TFH"),list(font = 4),
     main = list("Euler Diagram Overlaps K4me3 V2"))

plot(f,
     quantities = TRUE,
     fills = list(fill = c("#FF75A0", "#FCE38A", "#EAFFD0", "#95E1D3")),
     labels = c("TET","TFH"),list(font = 4),
     main = list("Euler Diagram Overlaps K4me3 V3"))
dev.off()
