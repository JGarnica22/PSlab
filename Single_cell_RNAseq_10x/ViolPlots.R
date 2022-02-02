## do graphs only with avglog2FC for every cluster

violplot <- function(pval,populations,comp){
# load all the excel file with the differential analysis of each populations and merge them
  for (i in populations){
    if (file.exists(paste0("out/DE_integrated_only_", i, "_KO.xlsx"))){
      x <- readxl::read_xlsx(paste0("out/DE_integrated_only_", i, "_KO.xlsx")) %>% 
        as.data.frame() %>% filter(p_val_adj < pval) %>% dplyr::select(gene, avg_log2FC)
      sc <- merge(sc, x, by = "gene", all = T)
      names(sc)[length(sc)] <- i
    }
  }
  # filter data and prepare for ggplot2
  sc[is.na(sc)] <- 0
  sc <- sc[rowSums(abs(sc[,-1]))>0,]
  sc <- tidyr::gather(sc, cluster, avglog2FC, -gene) %>% 
    filter(avglog2FC != 0)
  
  # do graph
  sc %>% ggplot(aes(avglog2FC, cluster)) +
    geom_jitter(size = 1)+
    geom_jitter(data = subset(sc, gene %in% features),
                color = "red")+
    geom_label_repel(data = subset(sc, gene %in% features == F),
                     aes(label = sc[sc$gene %in% features == F, "gene"]),
                     color = "navyblue",
                     size = 3,
                     force = 0.2,
                     # vjust = -0.1,
                     max.overlaps = 20)+
    geom_label_repel(data = subset(sc, gene %in% features),
                     aes(label = sc[sc$gene %in% features, "gene"]),
                     color = "red",
                     force = 0.25,
                     size = 3,
                     # vjust = -0.1,
                     max.overlaps = 20)+
    
    ylab("Cluster") +
    xlab("Avglog2FC DE genes") +
    theme(axis.ticks.y = element_text(size=15)
    )+
    titol(paste0("DE genes ", comp)) + theme_bw()
}  


