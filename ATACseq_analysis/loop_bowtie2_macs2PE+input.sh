#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J SANTAMARIA_16_loop
#BSUB -q bsc_ls
#BSUB -W 24:00
#BSUB -eo /gpfs/projects/cek26/SANTAMARIA_16/bsub_reports/SANTAMARIA_16_loop.err
#BSUB -oo /gpfs/projects/cek26/SANTAMARIA_16/bsub_reports/SANTAMARIA_16_loop.out
#BSUB -M 1800
#BSUB -n 64
#BSUB -R "span[ptile=16]"
#BSUB -x

# Set your working directory
wd=/gpfs/projects/cek26/SANTAMARIA_16
cd $wd

# Create a folder to store the fastq_files and upload all of them there
# mkdir fastq_files

# Create folders to store FastQC, trimmed sequences (if performing trimming) and alignment output
mkdir fastqc trimmed alignment alignment/no_dup to_bsub peak_calling

## process and get aligment of input files

module purge
module load java/1.8.0u66 intel/2017.4 impi/2017.4 MKL/2017.4 gcc/5.3.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip \
BOWTIE/2.4.2 SAMTOOLS/1.9

c=$(find ./fastq_files -name "Input1*2.fastq.gz" -exec basename {} \;)

java -jar ../software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 64 \
fastq_files/$(cut -d'_' -f1,2,3 <<< $c).fastq.gz fastq_files/$c  \
trimmed/trim_$(cut -d'_' -f1,2,3 <<< $c).fastq.gz trimmed/trim_1U_$(cut -d'_' -f1,2,3 <<< $c).fastq.gz \
trimmed/trim_$c trimmed/trim_2U_$c\
ILLUMINACLIP:../software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:8 MINLEN:15
bowtie2 -p 64 --very-sensitive -X 2000 -x ../genome_mus/chr_bowtie2/mm10 \
-1 trimmed/trim_$(cut -d'_' -f1,2,3 <<< $c).fastq.gz  -2 trimmed/trim_$c \
\| samtools view -@ 64 -u - \| samtools sort -@ 64 - \> alignment/bowtie_$(cut -d'_' -f1 <<< $c).bam
java -jar ../software/picard.jar MarkDuplicates -I alignment/bowtie_$(cut -d'_' -f1 <<< $c).bam -O alignment/no_dup/nd_$(cut -d'_' -f1 <<< $c).bam \
-M alignment/no_dup/$(cut -d'.' -f1 <<< $c)_log_dups.txt -REMOVE_DUPLICATES true
samtools index -@ 64 alignment/no_dup/nd_$(cut -d'_' -f1 <<< $c).bam

d=$(find ./fastq_files -name "Input2*2.fastq.gz" -exec basename {} \;)

java -jar ../software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 64 \
fastq_files/$(cut -d'_' -f1,2,3 <<< $d).fastq.gz fastq_files/$d  \
trimmed/trim_$(cut -d'_' -f1,2,3 <<< $d).fastq.gz trimmed/trim_1U_$(cut -d'_' -f1,2,3 <<< $d).fastq.gz \
trimmed/trim_$d trimmed/trim_2U_$d\
ILLUMINACLIP:../software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:8 MINLEN:15
bowtie2 -p 64 --very-sensitive -X 2000 -x ../genome_mus/chr_bowtie2/mm10 \
-1 trimmed/trim_$(cut -d'_' -f1,2,3 <<< $d).fastq.gz  -2 trimmed/trim_$d \
\| samtools view -@ 64 -u - \| samtools sort -@ 64 - \> alignment/bowtie_$(cut -d'_' -f1 <<< $d).bam
java -jar ../software/picard.jar MarkDuplicates -I alignment/bowtie_$(cut -d'_' -f1 <<< $d).bam -O alignment/no_dup/nd_$(cut -d'_' -f1 <<< $d).bam \
-M alignment/no_dup/$(cut -d'.' -f1 <<< $d)_log_dups.txt -REMOVE_DUPLICATES true
samtools index -@ 64 alignment/no_dup/nd_$(cut -d'_' -f1 <<< $d).bam


for f in $(find ./fastq_files -name "T*2.fastq.gz" -exec basename {} \;)
do
{
echo \#!/bin/bash
echo \#BSUB cwd $wd
echo \#BSUB -J atacseq_$(cut -d'_' -f1 <<< $f)
echo \#BSUB -q bsc_ls
echo \#BSUB -W 20:00
echo \#BSUB -eo $wd/bsub_reports/CHIP_$(cut -d'_' -f1 <<< $f).err
echo \#BSUB -oo $wd/bsub_reports/CHIP_$(cut -d'_' -f1 <<< $f).out
echo \#BSUB -M 1800
echo \#BSUB -n 64
echo \#BSUB -R "span[ptile=16]"
echo \#BSUB -x

echo module purge
echo module load java/1.8.0u66 intel/2017.4 impi/2017.4 MKL/2017.4 gcc/5.3.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip \
BOWTIE/2.4.2 SAMTOOLS/1.9
echo cd $wd
echo java -jar ../software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 64 \
fastq_files/$f trimmed/trim_$f \
echo java -jar ../software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 64 \
fastq_files/$(cut -d'_' -f1,2,3 <<< $f).fastq.gz trimmed/trim_$(cut -d'_' -f1,2,3 <<< $f).fastq.gz \
ILLUMINACLIP:../software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:8 MINLEN:15
echo bowtie2 -p 64 --very-sensitive -X 2000 -x ../genome_mus/chr_bowtie2/mm10 \
-1 trimmed/trim_$(cut -d'_' -f1,2,3 <<< $f).fastq.gz  -2 trimmed/trim_$f \
\| samtools view -@ 64 -u - \| samtools sort -@ 64 - \> alignment/bowtie_$(cut -d'_' -f1 <<< $f).bam
echo java -jar ../software/picard.jar MarkDuplicates -I alignment/bowtie_$(cut -d'_' -f1 <<< $f).bam -O alignment/no_dup/nd_$(cut -d'_' -f1 <<< $f).bam \
-M alignment/no_dup/$(cut -d'.' -f1 <<< $f)_log_dups.txt -REMOVE_DUPLICATES true
echo samtools index -@ 64 alignment/no_dup/nd_$(cut -d'_' -f1 <<< $f).bam
echo macs2 callpeak -t alignment/no_dup/nd_$(cut -d'_' -f1 <<< $f).bam -c alignment/no_dup/nd_$(cut -d'_' -f1 <<< $c).bam alignment/no_dup/nd_$(cut -d'_' -f1 <<< $d).bam \
 -f BAM -g 1.87e9 -q 0.05 --keep-dup all \
-n macs_$(cut -d'_' -f1 <<< $f) --outdir peak_calling 2\> peak_calling/macs_$(cut -d'_' -f1 <<< $f).log
} > to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
sed -i -e 's/\r$//' to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
bsub < to_bsub/atacseq_$(cut -d'_' -f1 <<< $f).sh
done

# Perform quality control FastQC tool for all processed samples and output a common pdf file
# Run FastQC for each .fastq.gz file and save results file in different folders in fastqc_analysis directory
# Load necessary modules
module purge
module load java/1.8.0u66 fastqc
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
