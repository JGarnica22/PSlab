# Single cell Multiome ATAC + Gene Expression analysis :dizzy:
In this tutorial we aim to describe the tools and steps needed to analyse Chromium Single Cell Multiome ATAC + Gene Expression sequencing data. This will allow us to study both gene expression and chromatin accessibility, and their relationship, since both measurements are taken on the very same cell.

As in other single-cell pipelines guides, there is a first part run in shell based on cellranger pipelines which basically involves alignment, filtering, barcode counting and peak calling. These should preferably be run in the cluster and a loop to run all cellranger steps for all the samples can be found in [cellranger_arc_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/cellranger_arc_loop.sh).
After this, this guideline recommends downloading cellranger outputs and importing them into Seurat package to do the downstream analysis and visualization. A script with Seurat functions to use is found at [10X_ATAC_RNA.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Multi-modal/10X_ATAC_RNA.R).

**Before starting** :heavy_exclamation_mark:: This guide assumes that you are already have the fastq files for analysis with `cellranger-arc`. If you are beginning with raw base call (BCL) files or more information about cellranger-arc visit [10X webpage](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc).  
</br>

## cellranger-arc count :rainbow:
**Guideline for cellranger-arc **v1.0***

### Get cellranger-arc-compatible reference
Before starting make sure to download the right cellranger-arc-compatible reference files (reference for human and mouse to donwload [here](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest)). You can also create a reference package using `cellranger-arc mkref` starting with a genome assembly FASTA file, a GTF file of gene annotations, and optionally a file of transcription factor motifs in JASPAR format.

### Create CSV file with libraries
Next, construct a 3-column CSV file that specifies the location of the ATAC and GEX FASTQ files associated with the sample. You can do this once you have uploaded your files into the cluster in the appropriate folder.

| Column name/parameter | Description                              |
|-----------------------|------------------------------------------|
| `fastqs`                | A fully qualified path to the *directory* containing the demultiplexed FASTQ files for this sample. This field does not accept comma-delimited paths. If you have multiple sets of fastq files for this library, add an additional row, and use the same `library_type` value. |
| `sample`                | Sample name. Important, sample name must have this format: SampleName_S1_L001_R1_001.fastq.gz |
| `library_type`          | This field is case-sensitive and must exactly match `Chromatin Accessibility` for a Multiome ATAC library and `Gene Expression` for a Multiome GEX library. |

For instance, a library CSV file would look like this:

````
fastqs,sample,library_type
/fastq_files/GEX_fastq,sample1,Gene Expression
/fastq_files/ATAC_fastq,sample1,Chromatin Accessibility
````  

### Run cellranger-arc count
NOTE: Run a separate instance of `cellranger-arc count` for each GEM well that was demultiplexed!

Run `cellranger-arc count` with the following arguments:

````
cellranger-arc count   --id=sample345 \
                       --reference=path/to/reference \
                       --libraries=arc_libraries.csv \
````

By default, `cellranger-arc` will use all the cores available on your system to execute pipeline stages. You can specify a different number of cores to use with the `--localcores` option. Similarly, `--localmem` will restrict the amount of memory (in GB) used.

The pipeline will create a new folder named with the **sample ID** you specified for its output. If this folder already exists, `cellranger-arc` will assume it is an existing pipestance and attempt to resume running it.

To see other command-line arguments run cellranger-arc count --help, here are some useful:
| Argument                      | Description                              |
|-------------------------------|------------------------------------------|
| `--gex-exclude-introns`       | Disable counting of intronic reads. In this mode we only count reads that are exonic and compatible with annotated splice junctions in the reference. Note: using this mode will reduce the UMI counts in the count matrix. |
| `--min-atac-count`            | **Cell-caller override**: define the minimum number of transposition events in peaks for a cell barcode. Note: this option must be specified in conjunction with `--min-gex-count`. If you specify `--min-atac-count=500` `--min-gex-count=300` then a barcode is considered a cell if it has at least 500 ATAC transposition events in peaks OR at least 300 GEX UMI counts. It is advisable to use these parameters only after reviewing the web summary generated using default parameters. |
| `--min-gex-count`             | **Cell-caller override**: define the minimum number of UMI counts for a cell barcode. Note: this option must be specified in conjunction with `--min-atac-count`. It is advisable to use these parameters only after reviewing the web summary generated using default parameters. |
| `--peaks`                     | **Peak-caller override**: specify peaks to use in downstream analyses from supplied BED file. Note that the file must only contain three columns specifying the contig, start, and end of the peaks with no comment lines. The peaks must not overlap each other. The file must be sorted by position with the same chromosome order as the reference package. |


### Output files
A successful cellranger-arc count run should conclude with a message similar to this:

````
Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!
 
yyyy-mm-dd hh:mm:ss Shutting down.
 
Saving pipestance info to "sample1/sample1.mri.tgz"
````

The output of the pipeline will be contained in the folder named with the sample ID you specified. The subfolder named `outs` will contain the main pipeline output files:

| File Name                     | Description                              |
|-------------------------------|------------------------------------------|
| `web_summary.html`            | Run summary metrics and charts in HTML format. |
| `summary.csv`                 | Run summary metrics in CSV format.       |
| `raw_feature_bc_matrix.h5`      | Raw feature barcode matrix stored as a [CSC sparse matrix](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/advanced/h5_matrices) in hdf5 format. The rows consist of all the gene and peak features concatenated together and the columns consist of all possible barcode sequences (numbering 736,320). |
| `raw_feature_bc_matrix`         |  Raw feature barcode matrix stored as a CSC sparse matrix in MEX format. The rows consist of all the gene and peak features concatenated together and the columns consist of all possible barcode sequences (numbering 736,320). |
| `per_barcode_metrics.csv`       | AATAC and GEX read count summaries generated for every barcode observed in the experiment. For more details see [Per-barcode metrics](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics). |
| `gex_possorted_bam.bam`         | GEX reads aligned to the genome and transcriptome annotated with barcode information in BAM format. |
| `gex_possorted_bam.bam.bai`     | Index for gex_possorted_bam.bam.         |
| `gex_molecule_info.h5`          | Count and barcode information for every GEX molecule observed in the experiment in hdf5 format. |
| `filtered_feature_bc_matrix.h5` |  Filtered feature barcode matrix stored as a CSC sparse matrix in hdf5 format. The rows consist of all the gene and peak features concatenated together (identical to raw feature barcode matrix) and the columns are restricted to those barcodes that are identified as cells. |
| `filtered_feature_bc_matrix`    | Filtered feature barcode matrix stored as a CSC sparse matrix in MEX format. The rows consist of all the gene and peak features concatenated together (identical to raw feature barcode matrix) and the columns are restricted to those barcodes that are identified as cells. |
| `cloupe.cloupe`                 | Loupe Browser visualization file with all the analysis outputs. |
| `atac_possorted_bam.bam.bai`    | ATAC reads aligned to the genome annotated with barcode information in BAM format. |
| `atac_possorted_bam.bam`        | Index for atac_possorted_bam.bam.        |
| `atac_peaks.bed`                | Locations of open-chromatin regions identified in this sample. These regions are referred to as "peaks". |
| `atac_peak_annotation.tsv`      | Annotations of peaks based on genomic proximity alone. Note that these are not functional annotations and they do not make use of linkage with GEX data. |
| `atac_fragments.tsv.gz`         | Count and barcode information for every ATAC fragment observed in the experiment in TSV format. |
| `atac_fragments.tsv.gz.tbi`     | Index for atac_fragments.tsv.gz.         |
| `atac_cut_sites.bigwig`         | Genome track of observed transposition sites in the experiment smoothed at a resolution of 400 bases in BIGWIG format. |
| `analysis`                      | Various secondary analyses that utilize the ATAC data, the GEX data, and their linkage: dimensionality reduction and clustering results for the ATAC and GEX data, differential expression, and differential accessibility for all clustering results above and linkage between ATAC and GEX data. See [Analysis Overview](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/analysis) for more information. |


## Aggregation of samples 
:question: cellranger-arc does not offer tools to aggregate data from differents replicates yet, this should be done in Seurat.  
</br>

## Seurat WNN analysis RNA+ATAC :twisted_rightwards_arrows:


