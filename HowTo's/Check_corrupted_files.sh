# This script is intented for checking if compressed files (.gz) are corrupted, for instance after cluster uploading.

# Integrity of the files will be found in the check_corrupted.out file.



#!/bin/bash
#BSUB -cwd /gpfs/projects/cek26
#BSUB -J check_corrupted
#BSUB -q sequential
#BSUB -W 1:00
#BSUB -eo /gpfs/projects/cek26/check_corrupted.err
#BSUB -oo /gpfs/projects/cek26/check_corrupted.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"



# Switch to your working directory
wd=/gpfs/projects/cek26
cd $wd

for f in $(find Calgary/Char-files -name  "*.gz")
do
gzip -t $f && echo ok || echo corrupted
done
