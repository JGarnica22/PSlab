# ATACseq analysis pipeline :icecream:

In this pipeline it will be described all the steps and tools needed to process ATACseq data, from fastq file to differential analysis and post-peak analysis.

## Quality control: FastQC
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

## Adapter removal: Trimmomatic
Currently, due to the ubiquitous use of Illumina’s Nextera
library for ATAC-seq, overrepresentation of Nextera
sequencing adapters is often observed and should be
removed for accurate read alignment, atlhough is optional and some sensibility may be lost. Low-quality bases can also be eliminated using these tools.

`Trimmomatic` performs a variety of useful trimming tasks for illumina paired-end and single ended data. As it runs in java you need to download the [binary](http://www.usadellab.org/cms/?page=trimmomatic) and then call with `java -jar` command. Next in the command you need to speciy if the task is pair-end (PE) or single-end(SE) indicate your input (fastq files) and output files and then trimming and options settings. For instance, a single end mode call would look like this:
````
java -jar <path to trimmomatic jar> SE <input> <output> <step 1> ...
````
Most trimming steps take one or more settings, delimited by ':' (a colon)

* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64

After trimming of you files is recommendable run again a FastQC analysis to check that trimming worked well.

## Alignment
After the assessment of sequence quality and proper trimming and filtering, reads need to be mapped to a reference genome. Bowtie2 and BWA-MEM are memory efficient and fast. 


### Bowtie
#### Genome index
The first step for alignment is generate or download the genome index for bowtie2. You can download a complete package of already indexed (sequence, annotation, indexes for multiple aligners including BWA-MEM and bowtie2) from [Illumina's iGenomes site](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Alternatively, you can build a ndw une using `bowtie2-build` using <your_genome.fa> sequence as only input. Bowtie2 index must include a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2.

Then you can proceed to align your reads using:
````
bowtie2 --very-sensitive -X 2000 -x <path_to_index/base> -U <your_file> | samtools view -u - | samtools sort - > <your_file.bam>
````
Bowtie2 output is a SAM file, which contains alignment information for each input read. The SAM can be compressed to a binary format (BAM) and sorted with SAMtools. This is best accomplished by piping the output from Bowtie2 directly to samtools.

### BWA-MEM
BWA works very similar to bowtie2, again we need to download and indicate the indexes of your reference genome, found at iGenomes too, or make them using `bwa index`. Next, we can procced with alignment using `bwa mem` command. Here is an example, like bowtie the output format is sam file, if you want to get bam files the best is to pipe it directly to samtools.
````
bwa mem -t -M genome_mus/BWA/index/GRCm38 <your_file> \
| samtools sort -@32 -o <aligned_file.bam> -
````
**-M**: Mark shorter split hits as secondary (for `Picard` compatibility)

## Remove PCR duplicates
PCR duplicates are exact copies of DNA fragments that arise during PCR. Since they are artifacts of the library preparation procedure, they may interfere with the biological signal of interest. One commonly used program for removing PCR duplicates is Picard’s `MarkDuplicates`.

````
java -jar path_to/picard.jar MarkDuplicates I=aligned_file.bam O=aligned_file_nodup.bam \
M=log_dups.txt REMOVE_DUPLICATES=true
````
M= allows to indicate a log file with duplication metrics
REMOVE_DUPLICATES= if TRUE duplicates are removed from output file, FALSE is the default option.

## PEAK CALLING
Peak calling is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of ATAC-seq experiment. There are various tools that are available for peak calling. One of the more commonly used peak callers is **MACS2** using the function `macs2 callpeak`, as in this example:
````
macs2 callpeak -t <file_to_call> -f <file_format> -q 0.05 --nomodel --extsize 150 --keep-dup all \
-n <prefix_output> --outdir <output_directory> 2\> <log_file>.log
````
Some details about options used:

**-t**: This is the only REQUIRED parameter for MACS. If you have more than one alignment file, you can specify them as -t A B C. MACS will pool up all these files together.
**-n**: The name string of the experiment. MACS will use this string NAME to create output files like NAME_peaks.xls, NAME_negative_peaks.xls, NAME_peaks.bed...
**--outdir**: MACS2 will save all output files into the specified folder for this option.
**-f**: Format of tag file can be ELAND, BED, ELANDMULTI, ELANDEXPORT, SAM, BAM, BOWTIE, BAMPE, or BEDPE. 
**-q**: The q-value (minimum FDR) cutoff to call significant regions. Default is 0.05.
**--nomodel**: to bypass building the shifting model.
**--extsize**: While --nomodel is set, MACS uses this parameter to extend reads in 5'->3' direction to fix-sized fragments.
**--keep-dup**: It controls the MACS behavior towards duplicate tags at the exact same location. Need to be set at `all` if duplicated were previously removed.

## Cluster loop :curly_loop:
Use loop scripts to perform all the steps aforementioned automatically and in parellel in the cluster. All folders and subfolders will be created, just make sure to have your fastq files in a fasq_files folder, and to add and indicate the directory of your reference genome and software (Trimmomatic and Picard).
