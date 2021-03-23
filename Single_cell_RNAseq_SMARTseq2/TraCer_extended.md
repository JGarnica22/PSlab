# TraCer :dog:

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

#### Assemble output

For each cell, an /<output_directory>/<cell_name> directory will be created. This will contain the following subdirectories.

1. 1. <output_directory>/<cell_name>/aligned_reads: This contains the output from Bowtie2 with the sequences of the reads that aligned to the synthetic genomes.
2. <output_directory>/<cell_name>/Trinity_output: Contains fasta files for each locus where contigs could be assembled. Also two text files that log successful and unsuccessful assemblies.
3. <output_directory>/<cell_name>/IgBLAST_output: Files with the output from IgBLAST for the contigs from each locus.
4. <output_directory>/<cell_name>/unfiltered_TCR_seqs: Files describing the TCR sequences that were assembled prior to filtering by expression if necessary.
    * unfiltered_TCRs.txt : text file containing TCR details. Begins with count of productive/total rearrangements detected for each locus. Then details of each detected recombinant.
    * <cell_name>_TCRseqs.fa : fasta file containing full-length, reconstructed TCR sequences.
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
* --loci : The specific loci that you wish to summarise. TraCeR knows about alpha (A), beta (B), gamma (G) and delta (D) for mouse and human TCRs. Other locus names can be specified if you use build. By default, TraCeR will attempt to summarise alpha and beta sequences. To change this, pass a space-delimited list of locus names. For example, to only look for gamma/delta sequences use --loci G D. Default = A B
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
* clonotype_network_[with|without]_identifiers.<graph_format>: graphical representation of clonotype networks either with full recombinant identifiers or just lines indicating presence/absence of recombinants.
