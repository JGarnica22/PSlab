# change direction of methylation and K27me3
# as they should be inversely correlated with rna expression and other depositions
col <- grep("K27me3|MEth",names(hm2))
hm3 <- hm2 %>%
  filter(abs(scRNAseq.TR1vsTFH.1) > 1)
  # filter(row.names(.) %in% table)
hm3[,col] <- hm3[,col]*(-1)
matrix <- as.matrix(hm3)
dendro <- as.dendrogram(hclust(d=dist(x=matrix)))
dendro_plot <- ggdendro::ggdendrogram(dendro, rotate = T)+
                theme(axis.text.y = element_blank())
ord <- order.dendrogram(dendro)
#order based on the dendogram
hm4 <- hm2 %>%
  filter(abs(scRNAseq.TR1vsTFH.1) > 1) %>% 
  # filter(row.names(.) %in% table) %>% 
        rownames_to_column(var="gene")
hm4$gene <- factor(hm4$gene, levels = hm4$gene[ord], 
                               ordered = TRUE) 
hm4 <- hm4 %>%
  mutate_if(is.numeric, scale, center=F) %>% 
  gather(comparison, fold, -gene)

hmpl2 <- hm4 %>%
  # filter(gene %in% table) %>% 
  ggplot(aes(comparison, gene, fill=fold))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "grey98", high = "red")+
  theme(
        # axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 60,  hjust=1),
        legend.position = "top")
grid.newpage()
print(hmpl2, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.92))

svg("figs/heatmap_chips_meth_ord.svg", width= 5, height = 10)
grid.newpage()
print(hmpl2, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 1.0))
dev.off()

svg("figs/heatmap_chips_meth_ord2.svg", width= 5, height = 13)
hmpl2
dev.off()
