# TrackViewer :eyeglasses:
This is a package with web interface for drawing elegant interactive tracks or lollipop plots to facilitate integrated analysis of multi-omics data. You can visualize mapped reads along with annotation as track layers for different NGS datasets such as ChIP-seq, RNA-seq, miRNA-seq, DNA-seq, SNPs and methylation data.  
<br/>

## Prepare your input files:
Input files for trackViewer can be **bigWig (.bw)** or **BED (.bed)** files.

You can obtain bw files from bam files with `bamCoverage`(deeptools) in Terminal. For that follow these instructions:  
<br/>

### 0. Install tools:
You will need `samtools` and `deeptools`, both included in conda (to install them see [HowTo_setupTerminalWLS](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalWLS.md) or [HowTo_SetupTerminalMac](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalMac.md) if you are working on Windows or Mac, respectively).  
<br/>

### 1. Create your folder to work:

````
cd /path/to/your/working/directory
mkdir bam_to_bw & cd "$_"
mkdir bam_files+index bw_files
````

Download bam files or move them into bam_files+index directory  
<br/>

### 2. Create an index file (.bai) for each bam file in the direcotry using `samtools` and then convert bam files into bigwigs with`bamCoverage` and store them in bw_files folder.

````
for f in $(find . -name "*.bam" -exec basename {} \;)
do
echo "Indexing:"$f
samtools index $f
echo $f".bai index file created"
echo Converting $f to bw
bamCoverage -b $f -o bw_files/Coverage_$f.bw -v
echo Coverage_$f.bw file created
done
````

Move bw files to the trackViewer _data_ folder in your computer:
````
mv bw_files/*.bw /path/to/your/project/data/folder #e.g. /Users/$USER/projects/trackViewer/data
````

\* In a WLS, you have to move it to **a Windows folder** for R analysis:
````
mv bw_files/*.bw /mnt/path/to/your/project/data/folder #e.g. /Users/$USER/projects/trackViewer/data /mnt/...
````


Use trackViewer: script.R
