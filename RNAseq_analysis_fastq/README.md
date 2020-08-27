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


### Alignment
