# Smart-Seq2 single-cell sequencing analysis :sunglasses:

This pipeline describes the analysis of Smart-seq2 data obtained from single cells. SMART stands for Switching Mechanism At the end of the 5’-end of the RNA Transcript, a key property of the reverse transcriptase enzyme from Maloney Murine Leukemia Virus (MMLV). During reverse transcription this enzyme adds a few nucleotides, generally 2-5 cytosines, when it reaches the 5’ end of the RNA template. These extra nucleotides on the newly-synthesized cDNA act as a docking site for a complementary oligonucleotide that carries 3 riboguanosines at its 3’ end. The reverse transcriptase is then able to switch templates and synthesize the complementary cDNA strand. Overall, optimization of this technique has improved both the yield and the length of transcripts from single-cell cDNA libraries.

## Overview of the pipeline
This pipeline processes every cell separately during alignment and gene expression quantification (**STAR**) and then everything is pooled into a common matrix in which every row is a gene and every column a cell. This output can be later analyzed with **Seurat**. In parallel, T cell receptors (TCR) sequences of each cell are reconstructed using **TraCer** and linked to single-cell data in Seurat. Finally, quality control of fastq files is done using **FastQC**.

All these tools are used as loops to cover all the files, taking into account that in the case of STAR and TraCer pair-ended fastq files must be paired for analysis while for FastQC each fastq file is analyzed separately.

[STAR_TraCer_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/SmartSeq2/STAR_TraCer_loop.sh) script generate different jobs with respective loops to process many samples in short time, then the outcome is pooled into objects to be analyzed on Seurat.

### STAR :star:
This pipelines uses STAR in order to align fastq files and at the same time to quantify the levels of transcripts for every gene with `--quantMode GeneCounts`. If not previously created, you must download your genome and annotation files and generate a STAR index (see [RNAseq_analysis_fastq](https://github.com/patriciasolesanchez/PSlab/tree/master/RNAseq_analysis_fastq) pipeline).

Since there are at least 1 fastq file (2 if paired-ended) for each cell, the best option is to do a loop for every file (see [STAR_TraCer_loop.sh](https://github.com/patriciasolesanchez/PSlab/blob/master/SmartSeq2/STAR_TraCer_loop.sh)). When making the loop for pair-end analysis make sure that files called indeed correspond to the same sample.

### TraCer :dog:

TraCer reconstructs the sequences of rearranged and expressed T cell receptor genes from single-cell RNA-seq data. It then uses the TCR sequences to identify cells that have the same receptor sequences and so derive from the same original clonally-expanded cell. TraCer depends on Bowtie2, Trinity, IgBlast and Kallisto or Salmon.

In NordIII cluster all these dependencies are contained in a singularity image. Therefore, in this case to run tracer you must use this command:

````
singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif <Command> <Options>
````

For more details and installation instructions see [TraCer webpage](https://github.com/Teichlab/tracer#installation).

### Usage

Tracer has three modes: assemble, summarise and build. We will be indicating how to use assemble and summarise.

* **Assemble** takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs TCR sequences.
* **Summarise** takes a set of directories containing output from the assemble phase (each directory represents a single cell) and summarises TCR recovery rates as well as generating clonotype networks.
* **Build** creates new combinatorial recombinomes for species other than the inbuilt Human and Mouse.

#### TraCer: Assemble
Usage
````
tracer assemble [options] <file_1> [<file_2>] <cell_name> <output_directory>
````
**Main arguments:**

* <file_1> : fastq file containing #1 mates from paired-end sequencing or all reads from single-end sequencing.
* <file_2> : fastq file containing #2 mates from paired-end sequencing. Do not use if your data are from single-end sequencing.
* <cell_name> : name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels etc.
* <output_directory> : directory for output. Will be created if it doesn't exist. Cell-specific output will go into /<output_directory>/<cell_name>. This path should be the same for every cell that you want to summarise together.

**Other arguments:**

* -p/--ncores <int> : number of processor cores available. This is passed to Bowtie2, Trinity, and Kallisto or Salmon. Default=1.
* -s/--species : Species from which the T cells were derived. Options are Mmus or Hsap for mouse or human data. If you have defined new species using the build mode, you should specify the same name here. Default = Mmus.
* --receptor_name : If, for some reason, you've used a different receptor name when using build then you should also specify it here. Default = TCR
* --loci : The specific loci that you wish to assemble. TraCeR knows about alpha (A), beta (B), gamma (G) and delta (D) for mouse and human TCRs.  This argument must be followed by another command line option. For example: --loci G D --p 1.
* -r/--resume_with_existing_files : if this is set, TraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if TraCeR died partway through a step and you want to resume where it left off.
* --single_end : use this option if your data are single-end reads. If this option is set you must specify fragment length and fragment sd as below.
* --max_junc_len : Maximum permitted length of junction string in recombinant identifier. Used to filter out artefacts. May need to be longer for TCR Delta.
* --invariant_sequences: Custom invariant sequence file.




