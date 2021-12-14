## bigwigs split
library(trackViewer)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggpubr)


f <- "Il10" # test first for only one gene
bw <- list.files(path= "data/aggr_multiome/subset_gex",
                 pattern= "Cluster*") # define the list of bigwig files to be used

cols2 <- c("#00B0F6", "#00BF7D", "#F8766D", "#A3A500", "#E76BF3") # define pallete of colors
#DefaultAssay(scdata) <- "ATAC"
#Idents(scdata) <- "predicted.id"

#Bigwig function to represent
bwig <- function(bw,w,tec,dir){
  pl <-  BigwigTrack(
    region = agr[[1]],
    bigwig = paste0("data/aggr_multiome/", dir, bw[w]),
    smooth = 3000,
    type = "coverage",
    y_label = tec,
    max.downsample = 3000, 
    downsample.rate = 0.1) +
    ggtitle(strsplit(bw[w], ".", fixed=T)[1][[1]][1]) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title.y = element_text(size = 9)) +
    xlab("")
  pl$layers[[1]]$aes_params$fill <- cols[w]
  return(pl)
}
pdf("figs/coverage_tracks_3.pdf", width = 16, height = 25)
for (f in features[-22]) {
  id <- get(f, org.SYMBOL2EG,)
  gr <- genes(TxDb, single.strand.genes.only=FALSE)[id]
  agr <- resize(gr, width(gr)[[1]]+100000, fix = "center")
  

  # genesinrange <- mapRangesToIds(TxDb, agr, type = "gene", ignore.strand = T)
  # tracks <- sapply(genesinrange[[1]][[1]],
  #                  function(z) {
  #                    track <- geneTrack(z,TxDb)[[1]]
  #                    return(track)
  #                  })
  

  rn <- list()
  rn_ctl <- list()
  rn_t <- list()
  ac <- list()
  ac_ctl <- list()
  ac_t <- list()
  for (w in 1:length(bw)) {
    rn[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/")
    ac[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/")
    rn_ctl[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/control/" )
    ac_ctl[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/control/")
    rn_t[[w]] <-  bwig(bw,w,"RNAseq", "subset_gex/treated/")
    ac_t[[w]] <- bwig(bw,w,"ATACseq", "subset_atac/treated/")
  }
  axy <- lapply(rn, function(x){layer_scales(x)$y$range$range[2]})
  maxy_rna <- axy %>% unlist() %>% max() %>% ceiling()
  axy <- lapply(ac, function(x){layer_scales(x)$y$range$range[2]})
  maxy_ac <- axy %>% unlist() %>% max() %>% ceiling()

  ann <- AnnotationPlot(scdata,region=GRangesToString(agr))
  ann$layers[[4]]$aes_params$size <- 6
  for (w in 1:length(bw)){
    rn[[w]] <- rn[[w]] + ylim(0,maxy_rna)
    rn_ctl[[w]] <- rn_ctl[[w]] + ylim(0,maxy_rna)
    rn_t[[w]] <- rn_t[[w]] + ylim(0,maxy_rna)
    ac[[w]] <- ac[[w]] + ylim(0,maxy_ac)
    ac_ctl[[w]] <- ac_ctl[[w]] + ylim(0,maxy_ac)
    ac_t[[w]] <- ac_t[[w]] + ylim(0,maxy_ac)
  }

  plt <- ggarrange(rn[[1]] , ac[[1]], rn[[2]] , ac[[2]],
                   rn[[3]], ac[[3]], rn[[4]], ac[[4]],
                   rn[[5]], ac[[5]],
                   ann, ncol=1)
  plt <- annotate_figure(plt,
                  top = text_grob("CD1d and CTL", color = "grey",
                                  face = "bold", size = 16))
  plt_ctl <- ggarrange(rn_ctl[[1]], ac_ctl[[1]], rn_ctl[[2]], ac_ctl[[2]],
                   rn_ctl[[3]], ac_ctl[[3]], rn_ctl[[4]], ac_ctl[[4]],
                   # rn_ctl[[5]], ac_ctl[[5]],
                   plot_spacer() + theme_transparent(),
                   plot_spacer() + theme_transparent(),
                   ann, ncol=1)
  plt_ctl <- annotate_figure(plt_ctl,
                         top = text_grob("CTL", color = "grey",
                                         face = "bold", size = 16))
  plt_t <- ggarrange(rn_t[[1]], ac_t[[1]], rn_t[[2]], ac_t[[2]],
                       rn_t[[3]], ac_t[[3]], rn_t[[4]], ac_t[[4]],
                       rn_t[[5]], ac_t[[5]], ann, ncol=1)
  plt_t <- annotate_figure(plt_t,
                             top = text_grob("CD1d", color = "grey",
                                             face = "bold", size = 16))
plt_all <- ggarrange(plt_ctl, plt_t, plt, ncol = 3)
plt_all <- annotate_figure(plt_all,
                       top = text_grob(f, color = "black",
                                       face = "bold", size = 19)) 
print(plt_all)
}
dev.off()
