#!/bin/bash

# Load modules in order to run R
module purge
module load gcc/7.2.0 R/3.6.3

# Merge all ReadPerGene dataframe column 2 (unstranded) into a unique dataframe for analysis with R, once all jobs are done.
R < Merge_Reads.R --save

# Remove unneeded files
cd mouse
rm alignment_counts/*SJ.out.tab alignment_counts/*.out.mate1
rm -r alignment_counts/*_STARgenome alignment_counts/*_STARpass1

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
mkdir alignment_counts/$(cut -d'.' -f1 <<< $f)_Logs 
mv alignment_counts/$(cut -d'.' -f1 <<< $f)*Log*.out alignment_counts/$(cut -d'.' -f1 <<< $f)_Logs
done
