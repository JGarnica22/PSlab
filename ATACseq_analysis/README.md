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

## Aligment
