# TrackViewer :eyeglasses:
This is a package with web interface for drawing elegant interactive tracks or lollipop plots to facilitate integrated analysis of multi-omics data. You can visualize mapped reads along with annotation as track layers for different NGS datasets such as ChIP-seq, RNA-seq, miRNA-seq, DNA-seq, SNPs and methylation data.  
<br/>

## Prepare your input files:
Input files for trackViewer can be **bigWig (.bw)** or **BED (.bed)** files.  


### bigWig files:

You can obtain bw files from bam files with `bamCoverage`(deeptools) in Terminal. For that, follow these instructions:  
<br/>

#### 0. Install tools:
You will need `samtools` and `deeptools`, both included in conda (to install them see [HowTo_setupTerminalWLS](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalWLS.md) or [HowTo_SetupTerminalMac](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalMac.md) if you are working on Windows or Mac, respectively).  
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
Bed files are easy to prepare


## Visualize data with trackViewer:
Once you have your input files prepared, you can use `trackViewer` to make your plots. 



