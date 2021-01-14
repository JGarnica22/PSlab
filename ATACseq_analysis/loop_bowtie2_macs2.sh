#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J atacseq_loop
#BSUB -q sequential
#BSUB -W 16:00
#BSUB -eo /gpfs/projects/cek26/ATACseq/bsub_reports/atacseq_loop.err
#BSUB -oo /gpfs/projects/cek26/ATACseq/bsub_reports/atacseq_loop.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x

# Load necessary modules
module purge
module load java/1.8.0u66 fastqc

# Set your working directory
wd=/gpfs/projects/cek26/ATACseq
cd $wd

# Create a folder to store the fastq_files and upload all of them there
# mkdir fastq_files

# Create folders to store FastQC, trimmed sequences (if performing trimming) and alignment output
mkdir fastqc trimmed alignment alignment/no_dup to_bsub peak_calling

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
{
echo \#!/bin/bash
echo \#BSUB cwd $wd
echo \#BSUB -J atacseq_$(cut -d'_' -f1 <<< $f)
echo \#BSUB -q bsc_ls
echo \#BSUB -W 8:00
echo \#BSUB -eo $wd/bsub_reports/atacseq_$(cut -d'_' -f1 <<< $f).err
echo \#BSUB -oo $wd/bsub_reports/atacseq_$(cut -d'_' -f1 <<< $f).out
echo \#BSUB -M 1800
echo \#BSUB -n 64
echo \#BSUB "span[ptile=16]"
echo \#BSUB -x

echo module purge
echo module load java/1.8.0u66 intel/2017.4 impi/2017.4 MKL/2017.4 gcc/5.3.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip \
BOWTIE/2.4.2 SAMTOOLS/1.9
echo cd /gpfs/projects/cek26/ATACseq
echo java -jar ../software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 64 \
fastq_files/$f trimmed/trim_$f \
ILLUMINACLIP:../software/Trimmomatic-0.39/adapters/NexteraSE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:8 MINLEN:15
echo bowtie2 -p 32 --very-sensitive -X 2000 -x ../genome_mus/bowtie2/genome -U trimmed/trim_$f \
\| samtools view -u - \| samtools sort - \> alignment/bowtie_$(cut -d'.' -f1 <<< $f).bam
echo java -jar ../software/picard.jar MarkDuplicates -I alignment/bowtie_$(cut -d'.' -f1 <<< $f).bam -O alignment/no_dup/nd_$(cut -d'.' -f1 <<< $f).bam \
-M alignment/no_dup/$(cut -d'.' -f1 <<< $f)_log_dups.txt -REMOVE_DUPLICATES true
echo samtools index alignment/no_dup/nd_$(cut -d'.' -f1 <<< $f).bam
echo macs2 callpeak -t alignment/no_dup/nd_$(cut -d'.' -f1 <<< $f).bam -f BAM -g 1.87e9 -q 0.05 --nomodel --shift 4 --extsize 150 --keep-dup all \
-n macs_$(cut -d'.' -f1 <<< $f) --outdir peak_calling 2\> peak_calling/macs_$(cut -d'.' -f1 <<< $f).log
} > to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
sed -i -e 's/\r$//' to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
bsub < to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
done

# Perform quality control FastQC tool for all processed samples and output a common pdf file
# Run FastQC for each .fastq.gz file and save results file in different folders in fastqc_analysis directory
files=fastq_files
fastqc=fastqc
for f in $(find ./$files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir $fastqc/$(cut -d'.' -f1 <<< $f)  
fastqc -t 64 $files/$f --extract -o $fastqc/$(cut -d'.' -f1 <<< $f)
$fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/summary.txt | \
	montage txt:- $fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title $(cut -d'.' -f1 <<< $f) $fastqc/$(cut -d'.' -f1 <<< $f).png
rm $fastqc/$(cut -d'.' -f1 <<< $f)/*.zip
done
# Create a pdf containing all images for all samples as a report
convert $fastqc/*.png $fastqc/fastqc_summary.pdf

# You can run again FastQC on trimmed files, evaluate and then proceed to aligment
while [ $(find ./trimmed -type f -name "*.fastq.gz" -print | wc -l) != $(find ./fastq_files -type f -name "*.fastq.gz" -print | wc -l) ]
do
sleep 10000;
done

files=trimmed
fastqc=fastqc_trim
mkdir $fastqc
for f in $(find ./$files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir $fastqc/$(cut -d'.' -f1 <<< $f)  
fastqc -t 64 $files/$f --extract -o $fastqc/$(cut -d'.' -f1 <<< $f)
$fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/summary.txt | \
	montage txt:- $fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title $(cut -d'.' -f1 <<< $f) $fastqc/$(cut -d'.' -f1 <<< $f).png
rm $fastqc/$(cut -d'.' -f1 <<< $f)/*.zip
done
convert $fastqc/*.png $fastqc/fastqc_summary.pdf

# Remove used scripts
rm -r to_bsub
