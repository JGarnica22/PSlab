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



## Seurat analysis :globe_with_meridians:

