# RNAseq data analysis pipeline :dna:

RNA-sequencing (RNAseq) generates raw data in the _fastq_ format. Each sample you sequence will generate a compressed _fastq.gz_ file. These files will be used as input for the following analysis pipeline.

The RNAseq analysis pipeline includes scripts to perform a `FastQC` quality control of each fastq file (sample) followed by the alignment of those to a reference genome and counting the number of reads per gene using the `STAR` tool. Note that this pipeline is thought to be run in the Cluster due to its high RAM requirements (for `STAR` ideally 32 GB); however, this code could be run in your computer Terminal, you will need longer time and your computer might overheat a bit. These scripts are written in `bash` language.  
</br>

## Quality control with FastQC
FastQC checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. FastQC is a cross-platform application, written in java. In theory it should run on any platform which has a suitable java runtime environment. To run it in the Cluster you will need to load the java and fastqc modules (see full bash script).

Visit their [website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

### FastQC input
You can use BAM, SAM or fastq files (any variant) as input files. We will be using compressed fastq files (_your_file.fastq.gz_), as they will be decompressed and analyzed subsequently.

### FastQC output
`FastQC` analysis provides some summary graphs to quicky evalute the quality of your data. Results can be exported to an HTML based report or to a PDF file. In this pipeline all graphs generated for all libraries are included in a single PDF file.

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

### Code
To run the `FastQC` tool with compressed fastq files you just need to use the following code:

````
fastqc your_file.fastq.gz --extract -o output_directory
````
</br>

To converge all graphs into one PDF file use:

````
summary.txt | \
montage txt:- /Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title your_file.png
````
</br>

## Read alignment with STAR
The first thing you need for an STAR alignment is creating a STAR index for the genome you are going to align you sequences to. This will need to be done just once, the first time you use that genome. Once you have created the index files, you can use them for later alignments. 

### STAR indexing

#### STAR indexing Input
In order to generate the STAR index you will need two files: your genome of reference and its respective GTF annotation (also needed later for counting reads). Ideally, these files should be updated and can be downloaded from GENCODE or Ensembl. To avoid compatibility problems during indexing, make sure that both files come from the same database.

***IMPORTANT** Remember that you cannot download files when in the Cluster, you must download them in your computer and then upload them using SSH (Cyberduck for instance).*  

Once downloaded, you must decompress the files by:
````
gzip -d your_files
````
</br>

Also if you genome is is in _.2bit_ format you will need to convert it to _.fa_ (fasta). To do so, use the 'twobittofa' tool. If not installed in your computer, you can install it through 'Conda'.
````
twoBitToFa your_genome.2bit your_genome.fa
````

#### Code
Finally you can generate your genome index by using:

````
STAR --runMode genomeGenerate --genomeDir /your_directory --genomeFastaFiles genome_reference.fa --sjdbGTFfile your_annotation.gtf --sjdbOverhang 49 --runThreadN 12
````
Note that this is action will take a lot of RAM. With `--runThreadN` you can assign the number of nodes working, the more nodes the faster it will be performed.
***IMPORTANT** Never use all of your computer nodes to run a task.

#### STAR index output
All STAR index files will be stored in the directory you indicated after `--genomeDir` option. All the files will be needed later for the alignment. 
***IMPORTANT** The file system needs to have at least 100GB of disk space available for a typical mammalian genome.*  
</br>

### Alignment :train:
After STAR index files are generated, the RNAseq reads (sequences) are mapped to the genome.

#### Alignment input
Here we will need the sequences to align in FASTA or FASTQ files format as well as all the files from STAR index. Also, when you want the number of reads per gene want to be counted while mapping (with `--quantMode GeneCounts` option), the annotation file used in the STAR index generation is required.

#### Code
````
STAR --runMode alignReads --genomeDir your_STARindex --readFilesIn your_seq_files --readFilesCommand zcat --outFileNamePrefix your_seq_file --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --sjdbGTFfile annotation_file.gtf --twopassMode Basic --quantMode GeneCounts --runThreadN 8 --genomeSAsparseD 3 --genomeSAindexNbases 12
````

#### Output
STAR produces multiple output files. All files have standard name, however, you can change the file prefixes using `--outFileNamePrefix /path/to/output/dir/prefix`. By default, this parameter is ./, i.e. all output files are written in the current directory.

* **Log files**: there are different types of log files generated. They provide progress statistics and information about the run, useful for quality control and troubleshooting.
* **Alignments in BAM format**: in our case we obtain output sorted by coordinate as *Aligned.sortedByCoord.out.bam* file. Other files can be obtained using the option `--outSAMtype`.
* **Reads per gene files**: STAR outputs read counts per gene into _ReadsPerGene.out.tab file_ with 4 columns which correspond to different strandedness options:

    * **column 1**: gene ID
    * **column 2**: counts for unstranded RNAseq
    * **column 3**: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
    * **column 4**: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

    You must select the output according to the [strandedness of your data](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html?showComment=1522859906203#c7586510008202176571) (for instance, <ins>standard Illumina would be unstranded</ins>).
    
## Running multiple jobs simultaneously in the Cluster :two_women_holding_hands:
Both STARindexing and aligments can take quite time to be performed due to the amount of data processed, and this depends on the amount of RAM available to use. In the Cluster you can highly increase the RAM used in the process in comparison to a computer and in case you have different replicates or samples, run them in parellel to save time. `Paralel_alignment_and_FastQC.sh` is a bash script intended to:
1. Run in parellel jobs the aligments while counting reads per gene for each sample or replicate you may have, and then convert each resulting bam file into Bigwig. To do so it creates an indepedent job script for each FASTQ file which is later removed.
2. Run a loop to perform a FasQC each of your samples and compile all graphs into a unique PDF file.

To use it, upload all necessary files in the proper directories and send the script as a job: `bsub < Paralel_alignment_and_FastQC.sh`

### Differential expression analysis

Once finished the alignment before [differential expression analysis](https://github.com/patriciasolesanchez/PSlab/tree/master/DE_analysis_RNAseq) we ought to put together all the data for all our samples in a single data frame. To do so you can use `Merge_Reads_R.sh` script. This script do not need to be sent as a job since is not memory-demanding, you can run it in the Cluster interface just by typing `Merge_Reads_R.sh` or `bash Merge_Reads_R.sh`. This script will run a R script (Merge_Reads.R) which need to be in the same directory and properly organize the output files as well as erase unnecessary files.

`Merge_Reads.R` basically perform a loop to concatenate the unstranded reads per gene (column 2) of each **ReadsPerGene.out.tab** file into a common txt file **All_reads.txt**.

