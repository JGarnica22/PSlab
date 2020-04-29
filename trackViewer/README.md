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

### 2. Create an index file (.bai) for each bam file in the direcotry using samtools
````
(conda activate samtools)
````
````
cd bam_files+index
````
````
for f in *.bam
do
echo "Indexing:"$f
samtools index $f
echo $f".bai index file created"
done
````

### 3. Convert bam files into bigwig (bw) files using `bamCoverage`
````
(conda activate deeptools)
````
````
for i in *.bam
do
echo "Converting "$i" to bw"
bamCoverage -b $i -o Coverage_$i.bw -v
echo "Coverage_"$i".bw file created"
done
````

Move files to bw_files folder
````
mv *.bw ~/projects/bam_to_bw/bw_files
````

In a WLS, you may want to move it to windows folder`
````
mv *bw /mnt/...
````


Use trackViewer: script.R
