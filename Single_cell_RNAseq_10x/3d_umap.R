# Maybe some libraries missing for this piece of code
library(scales)
library(plotly)
library(Seurat)

#Just to obtain the palette from ggplot2
cols <- hue_pal()(5)

scdata_origi2 <- scdata_origi
Embeddings(object = scdata_origi2, reduction = "umap")
scdata_origi2 <- RunUMAP(scdata_origi2,
                            dims = 1:10,
                            n.components = 3L) 
plot.data <- FetchData(object = scdata_origi2, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "celltype", "sampletype"))

figs <- list()
for (i in unique(plot.data$sampletype)) {
  tmp <- plot.data[plot.data$sampletype == i,]
  figs[[i]] <- plot_ly(data = tmp,
          x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
          color = ~celltype, 
          type = "scatter3d", 
          mode = "markers", 
          colors = cols,
          marker = list(size = 5, width=2)) %>%
        layout(title = paste0('UMAP 3D: ', i), plot_bgcolor = "#e5ecf6")
  names(figs[i]) <- i
          
}
print(figs)
