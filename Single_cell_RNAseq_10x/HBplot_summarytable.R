pops <- unique(merged$celltype2)
pval <- 0.05
# a object is a list of dataframes with differential analysis
l <- list()
for (pop in pops) {
  df <- data.frame(matrix(nrow=length(kos),ncol=5))
  #Set up df names
  names(df) <- c("Comparison","UP","ig.genesU", "DOWN","ig.genesD")
  df$Comparison <- kos
  for (j in 1:nrow(df)){
    ko <- df[j, 1]
    sc <- data.frame(gene = rownames(merged[["RNA"]]))
    populations <- unique(merged$celltype2[merged$sampletype == ko])
    for (i in populations){
      #Read the differential expression output tables ONLY for celltypes with more than 50 cells.
      if (ncol(merged[,merged$celltype2 == i & merged$sampletype == ko]) > 50) {
          x <- readxl::read_xlsx(paste0("outs/ko_vs_wt_foxp3/DE_CTLWT_vs_", ko, "_", i,"_wilcox.xlsx")) %>% 
            as.data.frame() %>% 
            dplyr::rename(avglog2FC = avg_diff) %>% 
            filter(p_val_adj < pval) %>% dplyr::select(gene, avglog2FC)
          sc <- merge(sc, x, by = "gene", all = T)
          names(sc)[length(sc)] <- i
      }
    }
    sc[is.na(sc)] <- 0
    if (ncol(sc) > 2) {
    sc <- sc[rowSums(abs(sc[,-1]))>0,]
    }
    sc <- tidyr::gather(sc, cluster, avglog2FC, -gene) %>% 
      filter(avglog2FC != 0)
    sc$avglog2FC <- - sc$avglog2FC
    sc <- sc[sc$cluster == pop,]
    scup <- sc %>% filter(avglog2FC > 0.5)
    df[j,2] <- nrow(scup)
    df[j,3] <- paste(sort(factor(scup$gene[scup$gene %in% features],
                      levels = features, ordered=TRUE)),
                      collapse=", ")
    scdown <- sc %>% filter(avglog2FC < -0.5)
    df[j,4] <- nrow(scdown)
    df[j,5] <- paste(sort(factor(scdown$gene[scdown$gene %in% features],
                      levels = features, ordered=TRUE)),
                      collapse=", ")
  }
  l[[pop]] <- df
}

writexl::write_xlsx(l, "outs/hbplots_tables_11JAN_2023.xlsx")
