#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J Loop_fastqc
#BSUB -q sequential
#BSUB -W 12:00
#BSUB -eo /gpfs/projects/cek26/experiments/all_rna/bsub_reports/Loop_fastqc.err
#BSUB -oo /gpfs/projects/cek26/experiments/all_rna/bsub_reports/Loop_fastqc.out
#BSUB -M 1800

# Load modules in order to run FastQC and STAR
module purge
module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/5.3.0 gcc/4.9.1 openssl/1.1.1c python/3.7.4_pip STAR/2.7.5a

# Set you working directory
wd=/gpfs/projects/cek26/experiments/all_rna

# Create a folder to store 'fastq_files' and upload all of them there
# also create directory before running the script for the reports
# mkdir fastq_files bsub_reports

# Create folders to store FastQC and STAR alignment
cd $wd
fastqc=fastqc
alignment_counts=alignment_counts
to_bsub=to_bsub
bsub_reports=bsub_reports

mkdir $fastqc $alignment_counts $to_bsub

# indicate location of STAR indexes and annotation file
STAR_index=/gpfs/projects/cek26/data/genomes/genome_mus/STAR_M25/
annotation=/gpfs/projects/cek26/data/genomes/genome_mus/gencode.vM25.annotation.gtf

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
{
echo \#!/bin/bash
echo \#BSUB cwd $wd
echo \#BSUB -J alignment_counts_$(cut -d'_' -f2 <<< $f)
echo \#BSUB -q bsc_ls
echo \#BSUB -W 4:00
echo \#BSUB -eo $wd/$bsub_reports/alignment_counts_$(cut -d'_' -f2 <<< $f).err
echo \#BSUB -oo $wd/$bsub_reports/alignment_counts_$(cut -d'_' -f2 <<< $f).out
echo \#BSUB -M 1800
echo \#BSUB -n 32
echo \#BSUB "span[ptile=16]"
echo \#BSUB -x

echo module purge
echo module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/5.3.0 gcc/4.9.1 openssl/1.1.1c python/3.7.4_pip STAR/2.7.5a SAMTOOLS/1.9 BEDTOOLS/2.25.0

echo cd $wd
echo STAR --runMode alignReads --genomeDir $STAR_index --readFilesIn fastq_files/$f --readFilesCommand zcat --outFileNamePrefix $alignment_counts/$(cut -d'_' -f2 <<< $f) --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $annotation --twopassMode Basic --quantMode GeneCounts --runThreadN 32

echo samtools index $alignment_counts/$(cut -d'_' -f2 <<< $f)\*.bam
echo bamCoverage -b $alignment_counts/$(cut -d'_' -f2 <<< $f)\*.bam -o $alignment_counts/Cover_$(cut -d'_' -f2 <<< $f).bw
} > $to_bsub/alignment_reads_$(cut -d'_' -f2 <<< $f).sh
sed -i -e 's/\r$//' $to_bsub/alignment_reads_$(cut -d'_' -f2 <<< $f).sh
bsub < $to_bsub/alignment_reads_$(cut -d'_' -f2 <<< $f).sh
done

# Perform quality control of RNAseq using FastQC tool for all processed samples and output a common pdf file
# Run FastQC for each .fastq.gz file and save results file in different folders in fastqc_analysis directory

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir $fastqc/$(cut -d'.' -f1 <<< $f)  
fastqc $fastq_files/$f --extract -o $fastqc/$(cut -d'.' -f1 <<< $f)
$fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/summary.txt | \
	montage txt:- $fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title $(cut -d'.' -f1 <<< $f) $fastqc/$(cut -d'.' -f1 <<< $f).png
rm $fastqc/$(cut -d'.' -f1 <<< $f)/*.zip
done

# Create a pdf containing all images for all samples as a report
convert $fastqc/*.png $fastqc/fastqc_summary.pdf

# Do not start next steps until previous pipelines are over
while [ $(find $alignment_counts -type f -name "*.bw" -print | wc -l) != $(find fastq_files -type f -name "*.fastq.gz" -print | wc -l) ]
do
sleep 300;
done

# Merge all samples gene counts into one dataframe using R
# Load modules in order to run R
module purge
module load gcc/7.2.0 R/3.6.3

# Merge all ReadPerGene dataframe column 2 (unstranded) into a unique dataframe for analysis with R, once all jobs are done.
{
echo RPG_list \<\- list.files\(path=paste0\(getwd\(\),\"/$alignment_counts\"\), pattern= \"*ReadsPerGene*\"\)
echo for \(g in 1\:length\(RPG_list\)\)\{
echo   x \<\- read.table\(paste0\(getwd\(\), \"/$alignment_counts/\", RPG_list\[g\]\)\) 
echo  x1 \<\- x[\-\(1\:4\),c\(1,2\)]
echo  names\(x1\) \<\- c\(\"Ensembl_id\", strsplit\(as.character\(RPG_list[g]\), split=\"_R\", fixed=TRUE\)[[1]][1]\)
echo if \(g == 1\)\{ all_reads \<\- x1 \} else \{all_reads \<\- merge\(all_reads, x1, by=\"Ensembl_id\"\) \}\}

echo rownames\(all_reads\) \<\- all_reads\$Ensembl_id
echo all_reads \<\- all_reads[,\-1]
      
echo write.table\(all_reads, \
            file = paste0\(getwd\(\),\"/$alignment_counts/Reads_all_samples.txt\"\), \
            sep = \"\\t\", quote = F, dec = \".\", row.names = T, col.names = T\) 
} > Merge_reads.R
R < Merge_reads.R --save

# Remove unneeded files
rm $alignment_counts/*SJ.out.tab $alignment_counts/*.out.mate1
rm -r $alignment_counts/*_STARgenome $alignment_counts/*_STARpass1

for f in $(find fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir $alignment_counts/$(cut -d'_' -f2 <<< $f)_Logs 
mv $alignment_counts/$(cut -d'_' -f2 <<< $f)*Log*.out $alignment_counts/$(cut -d'_' -f2 <<< $f)_Logs
done


