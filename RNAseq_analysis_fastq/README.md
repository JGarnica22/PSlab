# Fastq RNAseq analysis :sparkles:

In these scripts it will be described how to perform `fastQC` quality control of fastq files of RNAseq libraries as well as an alignment to a reference genome and counting the number of reads per gene using `STAR` tool. Note that these protocol are thought to be processed in the Cluster due to their highly RAM requeriments (for `STAR` ideally 32 GB), however, this code could be use in computer terminal. These scripts are written in `bash` language.

## Quality control of RNA-seq
FastQC checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. FastQC is a cross-platform application, written in java. In theory it should run on any platform which has a suitable java runtime environment. To run in the Cluster you will need to load the java and fastqc modules (see full bash script).

Visit their [website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

#### FastQC input
You can use here BAM, SAM or FastQ files (any variant). In this case we will be using compressed fastq files (your_file.fastq.gz), as they will be decompressed and analyzed subsequently. BAM files will be generated later too.

#### FastQC output
`FastQC` analysis provide some summary graphs to quicky evalute the quality of your data. Results can be exported to an HTML based report or to a PDF file. In this pipeline all graphs generated for all libraries are included in a single PDF file.

The graphs included are:
* Basic Statistics
* Per base sequence quality
* Per tile sequence quality
* Per sequence quality scores
* Per base sequence content
* Per sequence GC content
* Per base N content
* Sequence Length Distribution
* Sequence Duplication Levels
* Overrepresented sequences
* Adapter Content

#### Code
To run the `FastQC` from compressed fastq files you just need to use the following. However, this QC analysis is also included in the main script.

````
fastqc your_file.fastq.gz --extract -o output_directory
````
To converge all graph into on PDF file use:

````
summary.txt | \
montage txt:- /Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title your_file.png
````

## Read aligment with STAR

The first thing needed for STAR alignment is creating a STARindex for the genome you are going to align you sequences to, if you do not have it already. 

#### STARindex Input
In order to generate the STARindex you will need two files: your genome of reference and its respective GTF annotation (also needed later for counting reads). Ideally, these files should be as much updated as possible and can be downloaded from GENCODE or Ensembl, however, make sure that both files you will be using come from the same database, otherwise there could be problems.

***Note**: Remember that you cannot download files in the Cluster, you must download and then upload them using SSH (Cyberduck for instance). 

Once downloaded you must decompress the files by:
````
gzip -d your_files
````
Also if you genome is is in `.2bit` format you will need to convert it to `.fa`. To do so, install twobittofa with conda if not previously installed and use (this cannot be done on the Cluster):
````
twoBitToFa your_genome.2bit your_genome.fa
````

#### Code
Finally you can generate your genome index by using:

````
STAR --runMode genomeGenerate --genomeDir /your_directory --genomeFastaFiles genome_reference.fa --sjdbGTFfile your_annotation.gtf --sjdbOverhang 49 --runThreadN 12
````
Note that this is action will take a lot of RAM. With `--runThreadN` you can assign the number of nodes working, the more nodes the faster will be performed, however, never use all of your computer nodes for this task.

#### STARindex output

All STARindex files will be stored in the directory you indicated after `--genomeDir` option. All the files will be needed later for the alignment. 
****Important note**: The file system needs to have at least 100GB of disk space available for a typical mammalian genome.


### Alignment :train:
After STARindex files are generated the RNA-seq reads (sequences) are mapped to the genome.

#### Aligment input
Here we will need the sequences to align in FASTA or FASTQ files format as well as all the files from STARindex. Also, when number of reads per gene want to counted while mapping with `--quantMode GeneCounts` annotation file used in the STARindex generation is required.

#### Code
````
STAR --runMode alignReads --genomeDir your_STARindex --readFilesIn your_seq_files --readFilesCommand zcat --outFileNamePrefix your_seq_file --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --sjdbGTFfile annotation_file.gtf --twopassMode Basic --quantMode GeneCounts --runThreadN 8 --genomeSAsparseD 3 --genomeSAindexNbases 12
````

#### Output
STAR produces multiple output files. All files have standard name, however, you can change the file prefixes using `--outFileNamePrefix /path/to/output/dir/prefix`. By default, this parameter is ./, i.e. all output files are written in the current directory.

* **Log files**: there are different types of log files generated. They provides progress statistics and information about the run useful for quality control and troubleshooting.
* **Aligments in BAM format**: in our case we obtain output sorted by coordinate as *Aligned.sortedByCoord.out.bam* file. Other files can obteined in the opion `--outSAMtype`.
* **Reads per gene files**: STAR outputs read counts per gene into
ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:

    STAR outputs read counts per gene into
**ReadsPerGene.out.tab** file with 4 columns which correspond to different strandedness options:
    * **column 1**: gene ID
    * **column 2**: counts for unstranded RNA-seq
    * **column 3**: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
    * **column 4**: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

    You must select the output according to the [strandedness of your data](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html?showComment=1522859906203#c7586510008202176571) (for instance standard illumina would be unstranded).
    
## Multiple jobs generation in the Cluster :two_women_holding_hands:
Both STARindexing and aligments can take quite time to be performed due to the amount of data processed, and this depends on the amount of RAM available to use. In the Cluster you can highly increase the RAM used in the process in comparison to a computer and in case you have different replicates or samples, run them in parellel to save time. `Paralel_alignment_and_FastQC.sh` is a bash script intended to:
1. Run in parellel jobs the aligments while counting reads per gene for each sample or replicate you may have, and then convert each resulting bam file into Bigwig. To do so it creates an indepedent job script for each FASTQ file which is later removed.
2. Run a loop to perform a FasQC each of your samples and compile all graphs into a unique PDF file.

To use it, upload all necessary files in the proper directories and send the script as a job: `bsub < Paralel_alignment_and_FastQC.sh`

### Differential expression analysis

Once finished the alignment before [differential expression analysis](https://github.com/patriciasolesanchez/PSlab/tree/master/DE_analysis_RNAseq) we ought to put together all the date for all our samples in a single data frame. To do so you can use `Merge_Reads_R.sh` script. This script do not need to be sent as a job since is not memory-demanding, you can run it in the Cluster interface just by typing `Merge_Reads_R.sh` or `bash Merge_Reads_R.sh`. This script will run a R script (Merge_Reads.R) which need to be in the same directory and properly organize the output files as well as erase unnecessary files.

`Merge_Reads.R` basically perform a loop to concatenate the unstranded reads per gene (column 2) of each **ReadsPerGene.out.tab** file into a common txt file **All_reads.txt**.

