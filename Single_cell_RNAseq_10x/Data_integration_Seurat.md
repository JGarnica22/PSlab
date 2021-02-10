# Data integration :busts_in_silhouette:

When analyzing multiple single cell datasets obtained from different projects, techniques, species or on different days, it is very likely to see a batch effect, i.e. differences due to experimental artifacts and not biologically.

To solve that when analyzing such type of data `Seurat` offers methods for integration. These identify **'anchors'** between datasets to represent pairwise correspondences between individual cells, that is then hypothesized to orginate from the same biological state. These ‘anchors’ are then used to harmonize the datasets, or transfer information from one dataset to another.

In this guide we will show how to use this integration methods to then follow with your analyzis without batch effects.

## Merge your data (optional)
If your data from different datasets are already merged in the same file you can omit this step. Otherwise, you can combine two or more seurat object using `merge`:

````
combined <- merge(4k, y = 8k, add.cell.ids = c("4K", "8K"), project = "12K")
````

This will merge the raw count matrices of two Seurat objects and creates a new Seurat object with the resulting combined raw count matrix. To easily tell which original object any particular cell came from, you can set the `add.cell.ids` parameter with an c(x, y) vector, which will prepend the given identifier to the beginning of each cell name. The original project ID will remain stored in object meta data under `orig.ident`. To merge more than two Seurat objects, simply pass a vector of multiple Seurat objects to the y parameter.

## Preprocess your data
To construct a reference, we will identify ‘**anchors**’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element. Remember that for splitting your dataset, you must inlclude a column to metadata (in this exaple "project") indicating to which dataset belongs each cell.
````
list <- SplitObject(object, split.by = "project")
````

Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each, using variance stabilizing transformation ("vst").

````
for (i in 1:length(list)) {
    list[[i]] <- NormalizeData(list[[i]], verbose = FALSE)
    list[[i]] <- FindVariableFeatures(list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}
````

## Integration
Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input. 

````
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
````

We then pass these anchors to the `IntegrateData` function, which returns a Seurat object.

The returned object will contain a new `Assay`, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.

````
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
````

After running `IntegrateData`, the Seurat object will contain a new `Assay` with the integrated expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth.

We can then use this new integrated matrix for downstream analysis and visualization. Use this command to set the integrated matrix as the default asssay for analysis. 
````
DefaultAssay(integrated) <- "integrated"
````
**Note**: at this point data have already been normalized but it still needs to be scaled.

**IMPORTANT:** This integrated data is intented for PCA, t-SNE, UMAP and cluster analysis. To **find markers** you need to use raw (not integrated or original) data while preserving cluster labeling, as integration precisely mask differences between populations. This is why it is treated sample comparison as a two-step process. First, we integrate datasets together to consistently identify common cell populations across datasets. Second, we perform DE on the original (unintegrated) data, to find gene expression differences for each population. 

You can do so by using `DefaultAssay(integrated) <- "RNA"`
