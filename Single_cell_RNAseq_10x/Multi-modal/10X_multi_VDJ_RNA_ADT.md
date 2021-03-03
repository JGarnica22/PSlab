# Analysis of V(D)J, Gene Expression and Feature Barcode single-cell data :art:
The 5' Chromium Single Cell Immune Profiling Solution with Feature Barcode technology enables simultaneous profiling of V(D)J repertoire, cell surface protein, antigen specificity and gene expression data.

## cellranger multi :octopus:

The cellranger multi pipeline enables the analysis of these multiple library types together. The advantage of using the multi pipeline (as opposed to using `cellranger vdj` and `cellranger count` separately) is that it enables more consistent cell calling between the V(D)J and gene expression data. **Important: this guidelines are based on cellranger v6.0.**

**Note**: These guide lines assume that you already have fastq files either recieved from the facility or produced with `cellranger mkfastq`.

For more information visit [10X genomics webpage](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi).

### Run cellranger multi

To simultaneously generate single-cell feature counts, V(D)J sequences and annotations for a single library, run `cellranger multi`.

````
cellranger multi --id=samples345 --csv=project/samples345.csv
````

* `--id`: A unique run ID string that is also the output folder name.
* `--csv`: Path to config CSV file enumerating input libraries and analysis parameters.

The multi config CSV contains both the library definitions and experiment configuration variables. It is composed of up to four sections: `[gene-expression]`, `[feature]`, `[vdj]`, and `[libraries]`. The `[gene-expression]`, `[feature]`, and `[vdj]` sections have at most two columns, and are responsible for configuring their respective portions of the experiment. The `[libraries]` section specifies where input FASTQ files may be found. Click here for a multi config [template](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/multi_config_template.csv) and [example](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/multi_config_example.csv).

Basically, in the multi-config.csv file it must indicated reference genome for gene expression and vdj respective analysis, location of the fastq files and the features included in the libraries (*‘Gene Expression’, ‘VDJ’, ‘VDJ-T’, ‘VDJ-B’, ‘Antibody Capture’, ‘CRISPR Guide Capture’, or ‘Antigen Capture’*).

Optional arguments can also be added, such as previously seen --force-cells. For more information visit [10X genomics](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi).


## cellranger multi output files

A successful cellranger multi run should conclude with a message similar to this:

````
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!
 
yyyy-mm-dd hh:mm:ss Shutting down.
Saving pipestance info to "samples365.mri.tgz"
````

The output of the pipeline will be contained in a folder named with the **run ID** you specified (e.g. samples345). The subfolder named outs will contain the main pipeline outputs with the following structure:

![](https://support.10xgenomics.com/img/single-cell-vdj/multi-output-dir-structure.png)

The files in the `multi` folder contain raw data, i.e., the data for all barcodes including cells and background, while the files in the `per_sample_outs` directory are the filtered files, i.e. for cells in this sample. The pipeline also outputs a copy of the input VDJ reference in the `vdj_reference` directory.

| File Name            | Description                              |
|----------------------|------------------------------------------|
| `web_summary.html` | Run summary metrics and charts in HTML format |
| `metrics_summary.csv`  | Run summary metrics in CSV format        |
| `config.csv`           | The input multi config CSV               |
| `count`                | The results of any gene-expression and feature barcode analysis, similar to cellranger count |
| `vdj_b`                | The results of any V(D)J Immune Profiling analysis for any B cells, similar to cellranger vdj |
| `vdj_t`                | The results of any V(D)J Immune Profiling analysis for any T cells, similar to cellranger vdj |



***Note**: The gene expression library is representative of the entire pool of poly-adenylated mRNA transcripts captured within each partition (droplet). The TCR or BCR transcripts are then selectively amplified to create the V(D)J library. Therefore, the gene expression library has more power to detect partitions containing cells compared to the V(D)J library. If the multi pipeline is run with both gene expression and VDJ data, then barcodes which are not called as cells by using the gene expression data will be deleted from the V(D)J cell set.


## Agregate data with cellranger aggr
Many experiments involve generating multiple 10x libraries processed through different Gel Bead-in Emulsion (GEM) Wells on the Chromium instrument. Depending on the experimental design, these could be replicates from the same set of cells, cells from different tissue/time points from the same individual, or cells from different individuals. The `cellranger count`, `cellranger vdj`, and `cellranger multi` pipelines process data from a single GEM well. The aggr pipeline aggregates the outputs from multiple runs of `cellranger count/vdj/multi` and performs analysis on the combined data.
`cellranger aggr` is not designed for combining multiple sequencing runs from the same GEM Well (i.e., resequenced libraries). For that, pass the FASTQ files from multiple sequencing runs of the same GEM well to the count, vdj, or multi pipeline, as appropriate.

### Aggregating outputs from cellranger multi
The `cellranger aggr` command can take a CSV file specifying a list of `cellranger multi` output directories, and perform aggregation on any combination of 5' Gene Expression, Feature Barcode (cell surface protein, antigen, or CRISPR) and V(D)J libraries that are present in the individual runs of `cellranger multi`.

You can run the aggr pipeline as follows:
````
cellranger aggr --id=samples365_aggr --csv=aggr.csv
````

Create a CSV containing the following columns:

| Parameter   | Description                              |
|-------------|------------------------------------------|
| `sample_id`   | Unique identifier for this input GEM well. This will be used for labeling purposes only; it does not need to match any previous ID assigned to the GEM well. |
| `sample_outs` | Path to the per sample outs folder generated by cellranger multi. For example, if you processed your GEM well by calling cellranger multi --id=ID --csv=exp1.csv in some directory /DIR , and the sample was called Sample1, this path would be /DIR/ID/outs/per_sample_outs/Sample1 |
| `donor`       | An individual from whom adaptive immune cells (T cells, B cells) are collected (e.g. a sister and a brother would each be considered unique donors for the purposes of V(D)J aggregation). |
| `origin`      | The specific source from which a dataset of cells is derived. This could be a timepoint (pre- or post-treatment or vaccination or time A/B/C), a tissue (PBMC, tumor, lung), or other metadata (healthy, diseased, condition). Origins must be unique to each donor. Replicates (e.g. multiple libraries from the same population of cells) may share origins within a donor, which triggers additional replicate-based filtering. |

In addition to CSV columns described above, `cellranger aggr` accepts optional columns that may contain additional meta-data (e.g., vaccination status). These custom library annotations do not affect the analysis pipeline but can be visualized downstream in the Loupe V(D)J Browser.

`cellranger aggr` will auto-detect the presence of various libraries based on the structure and contents of the per sample outs folders. Apart from the change in the input CSV column (sample_outs instead of molecule_h5), the sections on aggregating outputs from cellranger count (depth normalization, batch correction etc.) applies here as well.

Here is an example of a config file:
````
library_id,library_outs,donor,origin,treatment
1-4,/media/josep/Elements/Parvus/cellranger/1-4/outs,1-4,1-4,1
1-5,/media/josep/Elements/Parvus/cellranger/1-5/outs,1-5,1-5,1
1-6,/media/josep/Elements/Parvus/cellranger/1-6/outs,1-6,1-6,1
1-7,/media/josep/Elements/Parvus/cellranger/1-7/outs,1-7,1-7,1
2-10,/media/josep/Elements/Parvus/cellranger/2-10/outs,2-10,2-10,2
2-11,/media/josep/Elements/Parvus/cellranger/2-11/outs,2-11,2-11,2
2-12,/media/josep/Elements/Parvus/cellranger/2-12/outs,2-12,2-12,2
2-6,/media/josep/Elements/Parvus/cellranger/2-6/outs,2-6,2-6,2
````

### cellranger aggr output files
Output files of `cellranger aggr` pipeline will be the same directories and files as individual jobs from `cellranger multi` but will be alocated in a single folder named as indiated in `--id`. `cellranger aggr` does not perform a cell-calling step, it simply aggregates the cell calls from each input job into a final set of cell calls.

## Seurat analysis :globe_with_meridians:

To analysis cellranger pipelines data on Seurat download `cellranger aggr` files, in case you have multiple GEM samples, otherwise download files directly from `cellranger multi` output. We recommend download and work with filtered matrixs.

Details on how to import and integrate different data type are detailed in [10X_RNA_VDJ_ADT.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/10X_RNA_ADT_VDJ.R)
