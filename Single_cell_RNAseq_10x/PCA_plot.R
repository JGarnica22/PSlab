library(tidyverse)
library(plyr)


## PCA
# `v` would be a list of dataframes with values to compare (either counts or log2FC...)
# skip this apply if NOT using DE data
la <- lapply(v, function(x){
                x <- x %>% filter(p_val_adj < 0.05) # exclude non-significant DE genes if your are using DE data
                names(x)[3] <- strsplit(x$comparison[1]," vs")[[1]][1] # change names of variables accordingly
                x <- x[,c(1,3)] # select only gene and log2FC column in case you are using DE data
                })
                
for (lg in c("all", "subset")){
le <- plyr::join_all(la, by="gene", type="full")
#if you want to perform for a subset of genes, or all...
if(lg == "subset"){
  le <- le %>% filter(gene %in% features)
}
le <- le %>% column_to_rownames(var="gene")
le[is.na(le)] <- 0 # needed when using log2FC only

pc <- prcomp(t(le))
pc_sum <- summary(pc)
PC1_varexpl <- pc_sum$importance[2,"PC1"]
PC2_varexpl <- pc_sum$importance[2,"PC2"]

# here indicate as column the rows depiciting the different categories to compare
plot <- as.data.frame(pc$x[,1:2]) %>% rownames_to_column(var="Var1") 

gp <- ggplot(plot, aes(x = PC1, y= PC2, 
                 color = Var1, label=Var1)) +
  geom_point(size= 5.5) +
  ggrepel::geom_text_repel(color="black", size=6)+
  ggprism::theme_prism()+
  labs(x= paste0("PC1 (", round(PC1_varexpl*100,1), "%)"),
       y= paste0("PC2 (", round(PC2_varexpl*100,1), "%)")) +
  ggtitle("Insert title")+ # insert plot title (optional)
  theme(axis.title.x = element_text(size = 20 ),
        axis.text.x = element_text(size = 20 ),
        axis.text.y = element_text(size = 20 ),
        axis.title.y = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(size=20, face="bold")
  )
        
svg("figs/PCA_plot.svg", width=6, height = 6) # save as svg
gp
dev.off()
