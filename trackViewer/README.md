# TrackViewer :eyeglasses:
This is a package with web interface for drawing elegant interactive tracks or lollipop plots to facilitate integrated analysis of multi-omics data. You can visualize mapped reads along with annotation as track layers for different NGS datasets such as ChIP-seq, RNA-seq, miRNA-seq, DNA-seq, SNPs and methylation data.  
<br/>

_BEFORE YOU USE THE TOOL..._

## Prepare your input files:
Input files for trackViewer can be **bigWig (.bw)** or **BED (.bed)** files.  
<br/>

### bigWig files:

bigWig is the format you will use mainly for sequencing data (to draw histograms with sequencing reads). You can obtain bw files from bam files (alignment files) with `bamCoverage`(deeptools) in Terminal. For that, follow these instructions:  
<br/>

#### 0. Install tools:
You will need `samtools` and `deeptools`, both included in conda (to install them see [HowTo_SetupTerminalWLS](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalWLS.md) or [HowTo_SetupTerminalMac](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalMac.md) if you are working on Windows or Mac, respectively).  
<br/>

#### 1. Create your folder to work:

Define your working directory:
````
wd=/path/to/your/project/.../trackViewer
#In WLS, you need to specify a Windows folder:
wd=/mnt/unit/$USER/.../trackviewer
````

Create new directory to work:
````
cd $wd
mkdir bam_to_bw && cd "$_"
````

Download bam files or move them into bam_to_bw directory.  
<br/>

#### 2. Create an index file (.bai) for each bam file in the direcotry using `samtools` and then convert bam files into bigwigs with `bamCoverage` and store them in data folder to be used by `RStudio`.

````
for f in $(find . -name "*.bam" -exec basename {} \;)
do
echo "Indexing:"$f
samtools index $f
echo $f".bai index file created"
echo Converting $f to bw
bamCoverage -b $f -o ../data/Coverage_$f.bw -v
echo Coverage_$f.bw file created
done
````
<br/>

### BED files:
BED (Browser Extensible Data) format will be used for data indicating selected regions, such as location of active enhacers or location of differentially methylated regions (DMRs). 

A BED file has the following format:

- BED lines have 3 required fields and nine additional optional fields.
- The number of fields per line must be consistent throughout any single set of data in an annotation track.
- The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used.  
<br/>

The BED fields are:
1. **chrom** - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).  
2. **chromStart** - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.  
3. **chromEnd** - The ending position of the feature in the chromosome or scaffold.  
<br/>

 _(optional fields)_  
   
4. **name** - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.  
5. **score** - A score between 0 and 100.  
6. **strand** - Defines the strand. Either "." (=no strand) or "+" or "-".  
7. **thickStart** - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.  
8. **thickEnd** - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).  
9. **itemRgb** - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.  
10. **blockCount** - The number of blocks (exons) in the BED line.  
11. **blockSizes** - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.  
12. **blockStarts** - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.  

_**Example**_
````
chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7  127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
````  
<br/>

BED files are easy to prepare, you can create them from a _.txt_ by doing the following:

1. Modify your file to contain the desired columns (usally for `trackViewer`, columns 1-5 are used). If using WLS you must do this from your Ubuntu terminal using your preferred text editor, for instance ´nano´:
````
nano your_file.txt
````
2. Eliminate column headers!
3. Change file extension from _.txt_ to _.bed_. From terminal use `mv` (*move*) command to change the extension of the file:
````
mv your_file.txt your_file.bed
````
<br/>

## Visualize data with trackViewer: :lollipop:
Once you have your input files prepared, you can use `trackViewer` to make your plots. Here we explain how to make two types of plots: chromosome views and lolliplots.

**Chromosome views** allow you to represent graphically in the same plot the different genomic and transcriptomic data from the same samples. In this case, _Chromosome_views.R_ is designed to represent data from **RNA-seq, ChIP, ATAC-seq, methylation** and inferred **active enhancers** (see [Active_enhancers_analysis folder] (https://github.com/patriciasolesanchez/PSlab/tree/master/Active_enhancers_analysis)) for your genes of interest.  

**Lolliplots** are a type of graph which allow you to clearly represent nucleotid-specific characteristics, such as methylation status or SNPs (Single Nucleotide Polymorphisms). _Lolliplots_overlapping.R_ script is prepared to represent the differential methylation between your samples in your genes of interest.  

Finally, _Chromosome_views_active_enhancers_guideslines_and_lolliplots.R_ script is thought to visualize the methylation status of a gene-associated active enhancers. It will make a first Chromosome view plot (like in _Chromosome_views.R_), but adding guidelines marking the position of active enhancers. Next, it will draw a lolliplot graph for each active enhancer (plotting the methylation status at the single C level in these active enhancers found 100 kb close to your genes of interest).
