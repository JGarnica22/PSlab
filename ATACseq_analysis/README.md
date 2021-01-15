# ATACseq analysis pipeline :icecream:

This pipeline describes all the steps and tools needed to process ATACseq data, from fastq files to differential analysis and post-peak analysis.  
</br>

## Quality control: FastQC
FastQC checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. FastQC is a cross-platform application, written in java. In theory it should run on any platform which has a suitable java runtime environment. To run it in the Cluster you will need to load the java and fastqc modules (see full bash script).

Visit their [website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

### FastQC input
You can use BAM, SAM or fastq files (any variant) as input files. We will be using compressed fastq files (_your_file.fastq.gz_), as they will be decompressed and analyzed subsequently.

### FastQC output
`FastQC` analysis provides some summary graphs to quickly evalute the quality of your data. Results can be exported to an HTML based report or to a PDF file. In this pipeline all graphs generated for all libraries are included in a single PDF file.

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
Currently, due to the ubiquitous use of Illumina’s Nextera library for ATAC-seq, overrepresentation of Nextera sequencing adapters is often observed and should be
removed for accurate read alignment, although is optional and some sensibility may be lost. Low-quality bases can also be eliminated using these tools.

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
</br>

## Alignment
After the assessment of sequence quality and proper trimming and filtering, reads need to be mapped to a reference genome. Bowtie2 and BWA-MEM are memory efficient and fast. 


### Bowtie
#### Genome index
The first step for alignment is to generate or download the genome index for bowtie2. You can download a complete package of already indexed (sequence, annotation, indexes for multiple aligners including BWA-MEM and bowtie2) from [Illumina's iGenomes site](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Alternatively, you can build a ndw une using `bowtie2-build` using <your_genome.fa> sequence as only input. Bowtie2 index must include a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2.

Then you can proceed to align your reads using:
````
bowtie2 --very-sensitive -X 2000 -x <path_to_index/base> -U <your_file> | samtools view -u - | samtools sort - > <your_file.bam>
````
Bowtie2 output is a SAM file, which contains alignment information for each input read. The SAM can be compressed to a binary format (BAM) and sorted with SAMtools. This is best accomplished by piping the output from Bowtie2 directly to samtools.

### BWA-MEM
BWA works very similar to bowtie2, again we need to download and indicate the indexes of your reference genome, found at iGenomes too, or make them using `bwa index`. Next, we can procced with the alignment using `bwa mem` command. Like bowtie the output format is a SAM file, if you want to get BAM files, the best is to pipe it directly to samtools. Here is an example:
````
bwa mem -t -M genome_mus/BWA/index/GRCm38 <your_file> \
| samtools sort -@32 -o <aligned_file.bam> -
````
**-M**: Mark shorter split hits as secondary (for `Picard` compatibility)  
</br>

## Remove PCR duplicates
PCR duplicates are exact copies of DNA fragments that arise during PCR. Since they are artifacts of the library preparation procedure, they may interfere with the biological signal of interest. One commonly used program for removing PCR duplicates is Picard’s `MarkDuplicates`.

````
java -jar path_to/picard.jar MarkDuplicates I=aligned_file.bam O=aligned_file_nodup.bam \
M=log_dups.txt REMOVE_DUPLICATES=true
````
M= allows to indicate a log file with duplication metrics
REMOVE_DUPLICATES= if TRUE duplicates are removed from output file, FALSE is the default option.  
</br>

## PEAK CALLING
Peak calling is a computational method used to identify areas in the genome that have been enriched with aligned reads and creating a 'peak' of reads. These correspond to accessible chromatin regions in an ATACseq experiment (equivalent to histone-marked or TF-binding regions in ChIPseq). There are various tools that are available for peak calling. One of the most used peak callers is **MACS2** using the function `macs2 callpeak`, as in this example:
````
macs2 callpeak -t <file_to_call> -f <file_format> -q 0.05 --nomodel --extsize 150 --keep-dup all \
-n <prefix_output> --outdir <output_directory> 2> <log_file>.log
````
Some details about options used:

* **-t**: This is the only REQUIRED parameter for MACS. If you have more than one alignment file, you can specify them as -t A B C. MACS will pool up all these files together. 

ENTENC QUE NOMÉS ESPECIFIQUES MÉS D'UN FILE A L'OPCIÓ -t QUAN SÓN LA MATEIXA MOSTRA? POTSER ES PODRIA ESPECIFICAR  

* **-n**: The name string of the experiment. MACS will use this string NAME to create output files like NAME_peaks.xls, NAME_negative_peaks.xls, NAME_peaks.bed...
* **--outdir**: MACS2 will save all output files into the specified folder.
* **-f**: Format of tag file can be ELAND, BED, ELANDMULTI, ELANDEXPORT, SAM, BAM, BOWTIE, BAMPE, or BEDPE. Default is “AUTO” which will allow MACS to decide the format automatically.
* **-q**: The q-value (minimum FDR) cutoff to call significant regions. Default is 0.05.
* **--nomodel**: to bypass building the shifting model.
* **--extsize**: While --nomodel is set, MACS uses this parameter to extend reads in 5'->3' direction to fix-sized fragments.
* **--keep-dup**: It controls the MACS behavior towards duplicate tags at the exact same location. Need to be set at `all` if duplicated were previously removed (i.e. with Picard).  
</br>

## Cluster loop :curly_loop:
Use loop scripts to perform all the steps aforementioned automatically and in parellel in the cluster. You have a _loop_bowtie2_macs2.sh_ and a _loop_bwa_macs2.sh_ script that will perform FastQC, Adapter trimming, Alignment and Peak calling; as indicated in the script name, they either use Bowtie2 or BWA-MEM as aligners.  

All folders and subfolders will be created automatically in the process, just make sure to have your fastq files in a fastq_files folder and to add and indicate the directory of your reference genome and software (Trimmomatic and Picard).

The following data analysis comparing peaks between samples and visualizing the results shall be performed using `R` language, and preferably locally, not in the cluster.  
</br>

## DiffBind :sunrise_over_mountains:
DiffBind provides functions for processing DNA data enriched for genomic loci including ChIPseq and ATACseq. It is designed to work with aligned sequence reads identified by a peak caller. The tool is optimized to work with multiple peak sets simultaneously and to identify sites that are differentially bound between sample groups.

As input for DiffBind you will need the peak files by macs2 (.xls) as well as the indexed alignment files (.bam and .bai) generated before.

Generally, data processing with DiffBind involves these phases:

### Reading in peaksets
Usually peaksets are derived from peak callers. The easiest way to read in peaksets is using a comma-separated value sample sheet or creating a dataframe with one line for each peakset. 

With your sample sheet you can generate your DBA object, which measures how many peaks are in each peakset, as well as (indicated in the first line) the total number of unique peaks after merging the overlapping ones. Also, a correlation heatmap and PCA can be generated, which gives an initial clustering of the samples using the cross-correlations of each row of the binding matrix.
````
dbdata <- dba(sampleSheet=<sample_sheet>)
plot(dbdadta)
dba.plotPCA(dbdata, label="ID")
````

### Counting reads
Once a consensus peakset has been derived, DiffBind can use the supplied sequence read files to count how many reads overlap each interval for each
unique sample. The final result of counting is a binding affinity matrix containing a read count for each sample at every consensus binding site, whether or not it was identified as a peak in that sample. Thus, a new column is added called **FRiP**, which stands for Fraction of Reads in Peaks. This is the proportion of reads for that sample that overlap a peak in the consensus peakset, and can be used to indicate which samples show more enrichment overall. With this matrix, the samples can be re-clustered using affinity, rather than occupancy, data. Same heatmap and PCA as previously can be generated to analyse this time the affinity scores.
````
dbdata <- dba.count(dbdata)
````

### Differential binding affinity analysis
Differential binding affinity analysis identifies significant differentially bound sites between sample groups. This will assign a p-value and FDR to each candidate binding site indicating the confidence that they are differentially bound.

First of all, before running the differential analysis, we need to tell DiffBind which cell lines fall in which groups based on our sample sheet, then we can perform the analysis.
````
dbadata <- dba.contrast(dbdata, categories =)
dbdata <- dba.analyze(dbdata)
````
This shows how many sites are identified as significant differentially bound (DB) using the default threshold of FDR≤0.05.
Again, we can perform PCA and Heatmaps with the affinity scores of these differentially bound sites. However, these plots are not results in the sense that the analysis is selecting for sites that differ between the two conditions, and hence are expected to form clusters.

### Reporting
Reporting mechanism enables differentially bound sites to be extracted for further processing, such as annotation, motif, and pathway analyses.
````
dbreport <- dba.report(dbdata, th = , fold = )
````
This command returns a GRanges object, appropiate for downstream processing. You can filter your report based on FDR threeshold (`th`) and/or Fold Change (`fold`).  
</br>

# Annotation :name_badge:
The next steps once we know all the peaks found in our samples and the differentially bound peaks between conditions is to know where these regions fall on the genome and next to what genes.

Here are described two annotations tools: `ChIPseeker` (run in R, from Bioconductor) and Homer's `annotatePeaks` (run in Terminal). We recommend using ChIPseeker directly in R after DiffBind analysis.  

## ChIPseeker :eyeglasses:
ChIPseeker is a very useful tool to annotate peak data analysis and visualize results in ggplot graphs.

### Input
ChIPseeker works with GRanges objects, so we can use the output of DiffBind reports or convert from our .bed or .xlsx files. To convert bed file to GRanges you can use `readPeakFile()` and to convert .xlsx files (MACS2 output) you can use:
````
peak <- read.delim(<peak file>, comment.char = "#") %>% toGRanges()
````

### Annotation
Afterwards, you can annotate your peaks using `annotatePeak()`. 
````
# Make sure your chr are in "chr1" format by:
seqlevelsStyle(peak) <- "UCSC"
# Annotate
peakAnno <- annotatePeak(peak, tssRegion=<Region Range of TSS>, TxDb=<TxDb object>, annoDb=<annotation package>)
````
This creates a csAnno object which can be easily exported as a data frame. Next, different graphs can be generated with this information.
````
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno, vennpie=T)
````

### Functional enrichment analysis
`ChIPseeker` can also perform functional enrichment analysis to identify predominant biological themes among these genes by incorporating biological knowledge provided by biological ontologies. For instance, Gene Ontology (GO).
````
enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
````

#### Script
Use `DiffBind_ChIPseeker.R` script to run all these steps combining DiffBind and ChIPseeker, obtaining the annotated differentially bound peaks and a comprehensive group of graphs visualizing the peaks dataset characteristics from each sample and from differentially bound peaks, which will be exported in pdf files.  
</br>

## Homer annotatePeaks
See requirements and steps for installation at http://homer.ucsd.edu/homer/introduction/install.html

### Install Homer
````
#!/bin/bash
wget <link with latest version>
perl /USER/homer/configureHomer.pl -install
export PATH=$PATH:/home/USER/HOMER/.//bin/
````

### Run annotatePeaks
First, you need to generate a BED file with your peaks to be annotated. In this case, the differentially bound peaks.

BED files should have minimum 6 columns (separated by TABs, additional columns will be ignored).
* Column1: chromosome (in `chr1` format)
* Column2: starting position
* Column3: ending position
* Column4: Unique Peak ID
* Column5: not used
* Column6: Strand (+/- or 0/1, where 0="+", 1="-")

This can be generated from the `dba.report` output using this code:
````
report <- as.data.frame(dbdata.DB) %>%
mutate(Unique_ID=row.names(report)) %>%
select(c(1:3, 12)) %>% 
mutate(seqnames = sapply("chr", paste0, seqnames))
report[,5:6] <- NA
write.table (report, "out/diffbind_report.bed",
             sep = "\t", dec = ".", quote = F, row.names = F, col.names = F)
````

Afterwards run `annotatePeaks.pl` in your Terminal. You need to install your genome of use before that in case you do not have it already:
````
#!/bin/bash
# Check available genomes to download
perl /home/USER/ATACseq/HOMER/.//configureHomer.pl -list
# Install your genome
perl /home/USER/ATACseq/HOMER/.//configureHomer.pl -install mm10
# Run annotation
annotatePeaks.pl diffbind_report.bed mm10 -go <dir for GO analysis> > homer_anno2.txt
````
The first two arguments, the peak file and genome, are required, and must be the first two arguments. `annotatePeaks.pl` also offers annotation enrichment analysis by specifying `-go GO output directory`.

Finally, you can import `annotatePeak.pl` data again into R and merge it with differential analysis report.
````
anno <- read.table("data/homer_anno.txt",
                   header = T,
                   sep = "\t",
                   quote = "", 
                   dec = ".")
names(anno)[1] <- "PeakID"  
comp <- merge(mutate(report, PeakID=row.names(report)), anno, by = "PeakID", all = T)
comp <- comp[, -c(2:4, 6)]
````
