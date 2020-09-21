#!/bin/bash

#BSUB cwd /gpfs/projects/cek26
#BSUB -J Loop_fastqc_SA_generate
#BSUB -q bsc_debug
#BSUB -W 2:00
#BSUB -eo /gpfs/projects/cek26/mouse/reports/Loop_fastqc_SA_generate.err
#BSUB -oo /gpfs/projects/cek26/mouse/reports/Loop_fastqc_SA_generate.out
#BSUB -M 1800


# Load modules in order to run FastQC and STAR
module purge
module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/5.3.0 gcc/4.9.1 openssl/1.1.1c python/3.7.4_pip STAR/2.7.5a

# Create a folder to store the fastq_files and upload all of them there
# mkdir fastq_files

# Create folders to store FastQC and STAR alignment
cd /gpfs/projects/cek26/mouse
mkdir fastqc
mkdir alignment_counts
mkdir ../to_bsub

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
{
echo \#!/bin/bash
echo \#BSUB cwd /gpfs/projects/cek26
echo \#BSUB -J Aligment_counts_$(cut -d'.' -f1 <<< $f)
echo \#BSUB -q bsc_ls
echo \#BSUB -W 2:00
echo \#BSUB -eo /gpfs/projects/cek26/mouse/reports/aligment_counts_$(cut -d'.' -f1 <<< $f).err
echo \#BSUB -oo /gpfs/projects/cek26/mouse/reports/aligment_counts_$(cut -d'.' -f1 <<< $f).out
echo \#BSUB -M 1800

echo module purge
echo module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/5.3.0 gcc/4.9.1 openssl/1.1.1c python/3.7.4_pip STAR/2.7.5a SAMTOOLS/1.9 BEDTOOLS/2.25.0

echo cd /gpfs/projects/cek26/mouse
echo STAR --runMode alignReads --genomeDir STAR_index --readFilesIn fastq_files/$f --readFilesCommand zcat --outFileNamePrefix alignment_counts/$(cut -d'.' -f1 <<< $f)_ --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --sjdbGTFfile gencode.vM25.annotation.gtf --twopassMode Basic --quantMode GeneCounts --runThreadN 32 --genomeSAsparseD 3 --genomeSAindexNbases 12

echo samtools index alignment_counts/$(cut -d'.' -f1 <<< $f)*.bam
echo bamCoverage -b alignment_counts/$(cut -d'.' -f1 <<< $f)*.bam -o alignment_counts/Cover_$(cut -d'.' -f1 <<< $f).bw
} > ../to_bsub/aligment_reads_$(cut -d'.' -f1 <<< $f).sh
sed -i -e 's/\r$//' ../to_bsub/aligment_reads_$(cut -d'.' -f1 <<< $f).sh
bsub < ../to_bsub/aligment_reads_$(cut -d'.' -f1 <<< $f).sh
done

# Perform quality control of RNAseq using FastQC tool for all processed samples and output a common pdf file
# Run FastQC for each .fastq.gz file and save results file in different folders in fastqc_analysis directory

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir fastqc/$(cut -d'.' -f1 <<< $f)  
fastqc ./fastq_files/$f --extract -o fastqc/$(cut -d'.' -f1 <<< $f)
fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/summary.txt | \
	montage txt:- fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title $(cut -d'.' -f1 <<< $f) fastqc/$(cut -d'.' -f1 <<< $f).png
rm fastqc/$(cut -d'.' -f1 <<< $f)/*.zip
done

# Create a pdf containing all images for all samples as a report
convert fastqc/*.png fastqc/fastqc_summary.pdf

# remove used scripts
rm -r ../to_bsub
