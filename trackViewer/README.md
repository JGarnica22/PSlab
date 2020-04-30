# TrackViewer :eyeglasses:
This is a package with web interface for drawing elegant interactive tracks or lollipop plot to facilitate integrated analysis of multi-omics data. You can visualize mapped reads along with annotation as track layers for NGS dataset such as ChIP-seq, RNA-seq, miRNA-seq, DNA-seq, SNPs and methylation data.

## Prepare your input files:
Input files for trackviewr are **bigwig (.bw)** files

To obtain bw files from bam files in terminal, follow this instructions:

### 1. Install tools:
For this you need `samtools` and `deeptools`, both included in conda (to install them see [HowTo_setupTerminalWLS](https://github.com/patriciasolesanchez/PSlab/blob/master/HowTo's/HowTo_SetupTerminalWLS.md)).

Create new directories to work:
````
cd ~bam_to_bw
mkdir bam_files+index bw_files
````

Download bam files or move them into bam_files+index directory

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

In a WLS, you may want to move it to a windows folder for R analysis:
````
mv *.bw /mnt/...
````


Use trackViewer: script.R
