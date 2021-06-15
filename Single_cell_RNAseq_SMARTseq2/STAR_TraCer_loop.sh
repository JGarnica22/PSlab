#!/bin/bash

#BSUB cwd /gpfs/projects/cek26
#BSUB -J Smart_generate
#BSUB -q sequential
#BSUB -W 48:00
#BSUB -eo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_generate.err
#BSUB -oo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_generate.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"


wd=/gpfs/projects/cek26/SMARTseq2_BDC_INS
in=/gpfs/home/cek26/cek26527/mouse
cd $wd
mkdir to_bsub
alig_folder=alignment5
tracer=tracer2
mkdir $alig_folder $tracer fastqc

# Split fastq files in groups to parallelize the jobs, here we will be splitting the files in 3:
part=$(expr $(ls fastq_files/



.fastq.gz | wc -l) / 3)

for i in \$Group_1 \$Group_2 \$Group_3
do

{
echo \#!/bin/bash
echo \#BSUB cwd /gpfs/projects/cek26
echo \#BSUB -J Smart_STAR_$i
echo \#BSUB -q bsc_ls
echo \#BSUB -W 48:00
echo \#BSUB -eo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_STAR_$i.err
echo \#BSUB -oo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_STAR_$i.out
echo \#BSUB -M 3000
echo \#BSUB -n 128
echo \#BSUB -R "span[ptile=16]"
echo \#BSUB -x

echo module purge
echo module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/7.2.0 openssl/1.1.1c \
python/3.7.4_pip STAR/2.7.5a

echo cd $wd
echo all=\$\(find ./fastq_files -name "*_1.fastq.gz" -exec basename {} "\;" \)
echo Group_1\=\$\(echo \"\$all\" \| head \-n $part\)
echo Group_2\=\$\(echo \"\$all\" \| tail \-n $part\)
echo Group_3\=\$\(echo \"\$all\" \| tail \-n $(expr $part \* 2 + 2) \| head \-n $(expr $part + 2)\)

# Loop for STAR alignment and gene count
echo for f in $i 
echo do
echo if \[ \-f $alig_folder\/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)_ReadsPerGene.out.tab \]\; then
echo echo \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) "done"
echo else

# Perform alignment
echo echo doing STAR \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)
echo STAR --runMode alignReads --genomeDir $in/STAR_index --readFilesIn fastq_files/\$f fastq_files/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)_2*.gz --readFilesCommand zcat \
--outFileNamePrefix $alig_folder\/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) \
--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $in/gencode.vM25.annotation.gtf \
--twopassMode Basic --quantMode GeneCounts --runThreadN 64
echo fi
echo done
} > to_bsub/Smart_STAR_$i.sh
sed -i -e 's/\r$//' to_bsub/Smart_STAR_$i.sh
bsub < to_bsub/Smart_STAR_$i.sh

 # Loop for TraCer
{
echo \#!/bin/bash
echo \#BSUB cwd /gpfs/projects/cek26
echo \#BSUB -J Smart_tracer_$i
echo \#BSUB -q bsc_ls
echo \#BSUB -W 48:00
echo \#BSUB -eo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_tracer_$i.err
echo \#BSUB -oo /gpfs/projects/cek26/SMARTseq2_BDC_INS/reports/Smart_tracer_$i.out
echo \#BSUB -M 3000
echo \#BSUB -n 128
echo \#BSUB -R "span[ptile=16]"
echo \#BSUB -x

echo module purge
echo module load singularity

echo cd $wd
echo all=\$\(find ./fastq_files -name "*_1.fastq.gz" -exec basename {} "\;" \)
echo Group_1\=\$\(echo \"\$all\" \| head \-n $part\)
echo Group_2\=\$\(echo \"\$all\" \| tail \-n $part\)
echo Group_3\=\$\(echo \"\$all\" \| tail \-n $(expr $part \* 2 + 2) \| head \-n $(expr $part + 2)\)

# Loop for TraCer
echo for f in $i 
echo do
echo if \[ \-d $tracer\/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) \]\; then
echo echo \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) "done"
echo else
echo echo doing tracer \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)
echo singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif tracer assemble -p 32 -s Mmus --loci A B --max_junc_len 50 --resume_with_existing_files \
fastq_files/\$f fastq_files/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)_2*.gz \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) $tracer
echo fi
echo done
} > to_bsub/Smart_tracer_$i.sh
sed -i -e 's/\r$//' to_bsub/Smart_tracer_$i.sh
bsub < to_bsub/Smart_tracer_$i.sh

done
 

module purge
module load java/1.8.0u66 fastqc

# Perform quality control of RNAseq using FastQC tool for all processed samples and output a common pdf file
# Run FastQC for each .fastq.gz file and save results file in different folders in fastqc_analysis directory
for f in $(find ./fastq_files -name "*.gz" -exec basename {} \;)
do
mkdir fastqc/$(cut -d "." -f1 <<< $f)  
fastqc ./fastq_files/$f --extract -o fastqc/$(cut -d'.' -f1 <<< $f)
fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/summary.txt | \
	montage txt:- fastqc/$(cut -d'.' -f1 <<< $f)/$(cut -d'.' -f1 <<< $f)_fastqc/Images/*.png \
	-tile x3 -geometry +0.1+0.1 -title $(cut -d'.' -f1 <<< $f) fastqc/$(cut -d'.' -f1 <<< $f).png
rm fastqc/$(cut -d'.' -f1 <<< $f)/*.zip
done

# Create a pdf containing all images for all samples as a report
convert fastqc/*.png fastqc/fastqc_summary.pdf

# Merge read counts and summarize tracer output
while [ $(find ./$alig_folder -type f -name "*Read*" -print | wc -l) != $(find ./fastq_files -type f -name "*_1.fastq.gz" -print | wc -l) ]
do
sleep 300;
done

module purge
module load gcc/8.4.0 openmpi/4.0.2 MKL/2017.4 R/4.0.4 


# Merge all ReadPerGene dataframe column 2 (unstranded) into a unique dataframe for analysis with R, once all jobs are done.
{
echo RPG_list \<\- list.files\(path=paste0\(getwd\(\),\"/$alig_folder\"\), pattern= \"*ReadsPerGene*\"\)
echo for \(g in 1\:length\(RPG_list\)\)\{
echo   x \<\- read.table\(paste0\(getwd\(\), \"/$alig_folder/\", RPG_list\[g\]\)\) 
echo  x1 \<\- x[\-\(1\:4\),c\(1,2\)]
echo  names\(x1\) \<\- c\(\"Ensembl_id\", strsplit\(as.character\(RPG_list[g]\), split=\"Re\", fixed=TRUE\)[[1]][1]\)
echo if \(g == 1\)\{ all_reads \<\- x1 \} else \{all_reads \<\- merge\(all_reads, x1, by=\"Ensembl_id\"\) \}\}

echo rownames\(all_reads\) \<\- all_reads\$Ensembl_id
echo all_reads \<\- all_reads[,\-1]
      
echo write.table\(all_reads, \
            file = paste0\(getwd\(\),\"/$alig_folder/Reads_all_samples.txt\"\), \
            sep = \"\\t\", quote = F, dec = \".\", row.names = T, col.names = T\) 
} > Smart.R
R < Smart.R --save

# Summarise TraCer
while [ $(ls $tracer| wc -l) != $(find ./fastq_files -type f -name "*_1.fastq.gz" -print | wc -l) ]
do
sleep 300;
done

module purge
module load singularity

singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif tracer summarise -p 16 -s Mmus $tracer

