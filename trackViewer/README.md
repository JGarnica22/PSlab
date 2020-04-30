# TrackViewer :eyeglasses:
This is a package with web interface for drawing elegant interactive tracks or lollipop plot to facilitate integrated analysis of multi-omics data. You can visualize mapped reads along with annotation as track layers for NGS dataset such as ChIP-seq, RNA-seq, miRNA-seq, DNA-seq, SNPs and methylation data.

## Prepare your input files:
Input files for trackviewr are **bigwig (.bw)** files

To obtain bw files from bam files in terminal, follow this instructions:

### 1. Install tools:
For this you need `samtools` and `deeptools`, both included in conda (to install them see [HowTo_setupTerminalWLS](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalWLS.md)).

Define your working directory:
````
wd=path/your/directory/.../tracviewer
wd=/mnt/unit/USER.../trackviewer
````
Create new directories to work:
cd $wd
mkdir bam_to_bw && cd "$_"
````

Download bam files or move them into bam_to_bw directory

### 2. Create an index file (.bai) for each bam file in the direcotry using `samtools` and then convert bam files into bigwigs with`bamCoverage` and store them in data folder to be used by `RStudio`.

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

Use trackViewer: script.R
