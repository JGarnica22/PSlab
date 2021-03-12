# SMARTseq2 single-cell sequencing analysis :sunglasses:

This pipeline describes the analysis of SMARTseq2 data obtained from single cells. SMART stands for Switching Mechanism At the end of the 5’-end of the RNA Transcript, a key property of the reverse transcriptase enzyme from Maloney Murine Leukemia Virus (MMLV). During reverse transcription this enzyme adds a few nucleotides, generally 2-5 cytosines, when it reaches the 5’ end of the RNA template. These extra nucleotides on the newly-synthesized cDNA act as a docking site for a complementary oligonucleotide that carries 3 riboguanosines at its 3’ end. The reverse transcriptase is then able to switch templates and synthesize the complementary cDNA strand. Overall, optimization of this technique has improved both the yield and the length of transcripts from single-cell cDNA libraries.

## Overview of the pipeline
This pipeline processes every cell separately during alignment and gene expression quantification (**STAR**) and then everything is pooled into a common matrix in which every row is a gene and every column a cell. This output can be later analyzed with **Seurat**. In parallel, T cell receptor (TCR) sequences of each cell are reconstructed using **TraCer** and linked to single-cell data in Seurat. Finally, quality control of fastq files is done using **FastQC**.

All these tools are run within a loop to cover all the files, taking into account that in the case of STAR and TraCer pair-end fastq files must be paired for analysis while for FastQC each fastq file is analyzed separately.

[STAR_TraCer_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/SmartSeq2/STAR_TraCer_loop.sh) script generates different jobs with respective loops to process many samples in short time, then the outcome is pooled into objects to be analyzed on Seurat.  
</br>

### STAR :star:
This pipelines uses the STAR aligner in order to both align fastq files and quantify the levels of transcripts for every gene with `--quantMode GeneCounts`. Unless previously created, you must download your genome and annotation files and generate a STAR index (see [RNAseq_analysis_fastq](https://github.com/patriciasolesanchez/PSlab/tree/master/RNAseq_analysis_fastq) pipeline for more information).

Since there are at least 1 fastq file (2 if paired-end) for each cell, the best option is to do a loop for every file (see [STAR_TraCer_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/SmartSeq2/STAR_TraCer_loop.sh)). When making the loop for pair-end analysis make sure that files called indeed correspond to the same sample.  
</br>

### TraCer :dog:

TraCer reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. TraCer depends on Bowtie2, Trinity, IgBlast and Kallisto or Salmon.

In NordIII cluster all these dependencies are contained in a singularity image. Therefore, in this case to run tracer you must use this command:

````
singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif <Command> <Options>
````

For more details and installation instructions see [TraCer webpage](https://github.com/Teichlab/tracer#installation).  

### Usage

TraCer has three modes: assemble, summarise and build. We will be indicating how to use assemble and summarise.

* **Assemble** takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs TCR sequences.
* **Summarise** takes a set of directories containing output from the assemble phase (each directory represents a single cell) and summarises TCR recovery rates as well as generates clonotype networks.
* **Build** creates new combinatorial recombinome for species other than the inbuilt Human and Mouse.  

#### TraCer: Assemble
````
tracer assemble [options] <file_1> [<file_2>] <cell_name> <output_directory>
````
**Main arguments:**

* <file_1> : fastq file containing #1 mates from paired-end sequencing or all reads from single-end sequencing.
* <file_2> : fastq file containing #2 mates from paired-end sequencing. Do not use if your data are from single-end sequencing.
* <cell_name> : name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels, etc.
* <output_directory> : directory for output. Will be created if it doesn't exist. Cell-specific output will go into /<output_directory>/<cell_name>. This path should be the same for every cell that you want to summarise together.

**Other arguments:**

* -p/--ncores <int> : number of processor cores available. This is passed to Bowtie2, Trinity, and Kallisto or Salmon. Default=1.
* -s/--species : Species from which the T cells were derived. Options are Mmus or Hsap for mouse or human data. If you have defined new species using the build mode, you should specify the same name here. Default = Mmus.
* --receptor_name : If, for some reason, you've used a different receptor name when using build then you should also specify it here. Default = TCR
* --loci : The specific loci that you wish to assemble. TraCeR knows about alpha (A), beta (B), gamma (G) and delta (D) for mouse and human TCRs. This argument must be followed by another command line option. For example: --loci G D --p 1.
* -r/--resume_with_existing_files : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.
* --single_end : use this option if your data are single-end reads. If this option is set you must specify fragment length and fragment sd as below.
* --max_junc_len : Maximum permitted length of junction string in recombinant identifier. Used to filter out artefacts. May need to be longer for TCR Delta.
* --invariant_sequences: Custom invariant sequence file.  


#### Assemble output

For each cell, an /<output_directory>/<cell_name> directory will be created. This will contain the following subdirectories.

1. <output_directory>/<cell_name>/aligned_reads: This contains the output from Bowtie2 with the sequences of the reads that aligned to the synthetic genomes.
2. <output_directory>/<cell_name>/Trinity_output: Contains fasta files for each locus where contigs could be assembled. Also two text files that log successful and unsuccessful assemblies.
3. <output_directory>/<cell_name>/IgBLAST_output: Files with the output from IgBLAST for the contigs from each locus.
4. <output_directory>/<cell_name>/unfiltered_TCR_seqs: Files describing the TCR sequences that were assembled prior to filtering by expression if necessary.
    * unfiltered_TCRs.txt : text file containing TCR details. Begins with count of productive/total rearrangements detected for each locus. Then details of each detected recombinant.
    * <cell_name>\_TCRseqs.fa : fasta file containing full-length, reconstructed TCR sequences.
    * <cell_name>.pkl: Python pickle file containing the internal representation of the cell and its recombinants as used by TraCeR. This is used in the summarisation steps.
5. <output_directory>/<cell_name>/expression_quantification: Contains Kallisto/Salmon output with expression quantification of the entire transcriptome including the reconstructed TCRs. When option --small_index is used, this directory contains only the output of the quantification with the small index (built from reconstructed TCRs and only a subset of the base transcriptome; see above).
6. <output_directory>/<cell_name>/filtered_TCR_seqs: Contains the same files as the unfiltered directory above but these recombinants have been filtered so that only the two most highly expressed from each locus are retained. This resolves biologically implausible situtations where more than two recombinants are detected for a locus. This directory contains the final output with high-confidence TCR assignments.


#### TraCer: Summarise
````
tracer summarise [options] <input_dir>
````

**Main argument**
<input_dir> : directory containing subdirectories of each cell you want to summarise.

**Options**
* -u/--use_unfiltered : Set this flag to use unfiltered recombinants for summary and networks rather than the recombinants filtered by expression level.
* -s/--species : Species from which the T cells were derived. Options are Mmus or Hsap for mouse or human data. If you have defined new species using the build mode, you should specify the same name here. Default = Mmus.
* --receptor_name : If, for some reason, you've used a different receptor name when using build then you should also specify it here. Default = TCR
* --loci : The specific loci that you wish to summarise. TraCeR knows about alpha (A), beta (B), gamma (G) and delta (D) for mouse and human TCRs. Other locus names can be specified if you use build. By default, TraCeR will attempt to summarise alpha and beta sequences. To change this, pass a space-delimited list of locus names.
* -i/--keep_invariant : TraCeR attempts to identify invariant TCR cells by their characteristic TCRA gene segments. By default, these are removed before creation of clonotype networks. Setting this option retains the invariant cells in all stages.
* -g/--graph_format : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.
* --no_networks : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.


#### Summarise Output
Output is written to the newly created directory <input_dir>/filtered_TCR_summary or <input_dir>/unfiltered_TCR_summary depending on whether the --use_unfiltered option was set.

The following output files are generated:

* TCR_summary.txt: Summary statistics describing successful TCR reconstruction rates and the numbers of cells with 0, 1, 2 or more recombinants for each locus.
* recombinants.txt: List of TCR identifiers, lengths, productivities, and CDR3 sequences (nucleotide and amino acid) for each cell. Note: It's possible for non-productive rearrangements to still have a detectable CDR3 if a frameshift introduces a stop-codon after this region. In these cases, the CDR3 nucleotides and amino acids are still reported but are shown in lower-case.
* reconstructed_lengths_TCR[A|B].pdf and reconstructed_lengths_TCR[A|B].txt: Distribution plots (and text files with underlying data) showing the lengths of the VDJ regions from assembled TCR contigs. Longer contigs give higher-confidence segment assignments. Text files are only generated if at least one TCR is found for a locus. Plots are only generated if at least two TCRs are found for a locus.
* clonotype_sizes.pdf and clonotype_sizes.txt: Distribution of clonotype sizes as bar graph and text file.
* clonotype_network_[with|without]\_identifiers.<graph_format>: graphical representation of clonotype networks either with full recombinant identifiers or just lines indicating presence/absence of recombinants.  
</br>


## Seurat analysis :raised_hands:
Seurat analysis for SMARTseq2 is based on [Seurat.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Seurat.R) script of 10X pipeline with few modifications here detailed. After applying these modifications you can run the rest of the commands to analyze and generate graphs.

These changes are included in [Seurat_Smartseq2.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_SMARTseq2/Seurat_SmartSeq2.R) script with an analysis example. If you require to work with data coming from different projects, you can find info about single-cell integration in [Data_integration_Seurat.md](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Data_integration_Seurat.md).

#### Import data
In this case, instead of having the 3 files coming from the cellranger pipeline that can be imported with the function `Read10X()`, we need to load the txt file obtained previously from STAR quantification with all the information using `read.table()`.

Note that as long as your single-cell data is arranged as a dataframe (in which each rowname is a gene and each column represents a cell) you can convert it to a Seurat object using `CreateSeuratObject()`. 

#### Change gene nomenclature (optional)
Depending on the genome and annotation file used, your genes may be labelled with ensembl_id instead of the gene name. Convert them to gene name with `biomaRt` using the following code. It is easier to do this before creating the Seurat object.

````
species <- "mouse"
en <- sapply(strsplit(rownames(df), ".", fixed = T), "[", 1) %>% as.data.frame()
names(en) <- "ensembl_gene_id"
if (species == "mouse"){
  TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  org.SYMBOL2EG <- org.Mm.egSYMBOL2EG
  org.SYMBOL <- org.Mm.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
} else {
  TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  org.SYMBOL2EG <- org.Hs.egSYMBOL2EG
  org.SYMBOL <- org.Hs.egSYMBOL
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
}
BM <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
            mart = ensembl, verbose = T)
enbm <- merge(en, BM, by="ensembl_gene_id")
row.names(df) <- make.names(enbm$external_gene_name, unique = T)
````

**IMPORTANT**: Also check the nomenclature of mitchondrial genes to see if they are `mt-` or `mt.` prefixed when calculating the QC metrics.
````
if (species == "mouse") {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^mt.")
} else {
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT.")
}
````


#### Add metadata
Once you have created you Seurat object you can add as much information as you want using the function `AddMetaData(object, metadata, col.name = NULL)` or more easily adding columns to your seurat object metada: `seuratobject@metadata$newcol`.

At this point it is essential to which group, condition, sample type, etc. each of the cells corresponds to. There are many ways to do that, using a loop, with a dataframe containing the information... This will allow us to analyze the data based on these parameters. For instance:
````
scdata$sampletype <- colnames(scdata)
for (i in 1:length(sampletype)){
  scdata@meta.data[grep(paste0("-", i), colnames(scdata)), "sampletype"] <- sampletype[i]
}
````
or

````
for (i in rownames(scdata@meta.data)){
  scdata@meta.data[i , "sampletype"] <- sub(".*?_", "", samples[ i ,"SAMPLE NAME"])
}
````

Moreover, at this point you must add also the **TCR clonotype information** obtained from `TraCer`. In this case we recommend using `merge`. With this you will be able to visualize how the clonotypes distribute within the clusters and sample types or conditions. 

````
# import TraCer data
tracer <- read.csv(paste0(getwd(),"/data/cell_data.csv"))
# Optionally divide the information in different columns in order to be available to subfiltrate. 
tracer <- tracer %>% separate(A_productive, c("A_produc_TRV", "A_produc_CDR3", "A_produc_TRJ"), "_", remove = F)
tracer <- tracer %>% separate(B_productive, c("B_produc_TRV", "B_produc_CDR3", "B_produc_TRJ"), "_", remove = F)

# Create a column in metadata with the cellnames in order to be able to merge
scdata@meta.data$cell_name <- rownames(scdata@meta.data)
scdata@meta.data <- merge(scdata@meta.data, tracer, by = "cell_name", all = T )
rownames(scdata@meta.data) <- scdata@meta.data$cell_name

# Delete cell name column and cells that produced a TCR sequence but were discarded when seurat object was created
scdata@meta.data <- scdata@meta.data %>% dplyr::select(-cell_name) %>% filter(is.na(sampletype) == F)
````

After these steps you can peform an usual Seurat analysis as detailed in [Seurat.R](https://github.com/patriciasolesanchez/PSlab/blob/master/Single_cell_RNAseq_10x/Seurat.R).

