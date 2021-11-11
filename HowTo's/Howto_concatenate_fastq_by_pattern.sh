#!/bin/bash


# This script concatenate fastq files obtained in different runs based on their filename.

#### IF SAME NAME, IN DIFFERENT FOLDERS (1, 2, 3, last the merged one) #########
for f in $(find -name "*.fastq.gz" -exec basename {} \; | sort -u)
do
cat $1/$f $2/$f > $3/merged_$f
done

echo Concatenation performed


#########################################################################################
#### IF NOT SAME NAME BUT SAME PREFIX ############

# Concatenate files
for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo $(cut -d "_" -f1 <<< $f)\*$(cut -d "_" -f3 <<< $f)  >> samples.txt
done

cd fastq_files
for f in $(sort -u ../samples.txt)
do
cat *$f*1.fastq.gz > merged/$f"_1.fastq.gz"
cat *$f*2.fastq.gz > merged/$f"_2.fastq.gz"
done

echo Concatenation performed
