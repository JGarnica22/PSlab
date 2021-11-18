#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J bamcoverage
#BSUB -q bsc_debug
#BSUB -W 4:00
#BSUB -eo /gpfs/projects/cek26/experiments/SANTAMARIA_31/bsub_reports/bamcoverage.err
#BSUB -oo /gpfs/projects/cek26/experiments/SANTAMARIA_31/bsub_reports/bamcoverage.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x


###############################################################################################
# THIS SCRIPT GENERATE .bw FILES OUT OF .bam FILES USING THE TOOLS bamCoverage from DEEPTOOLS
# .bw FILES CAN BE THEN USED FOR GENOMIC VIEWERS UTILITIES
################################################################################################

module purge
module load intel/2021.4 mkl/2021.4 python/3.7.12 # Careful with old and new Nord III cluster

# Set your working directory
wd=/gpfs/projects/cek26/experiments/SANTAMARIA_31/alignment/no_dup
cd $wd
mkdir bw

## indexing already done (.bai files present)
mkdir bw
for f in $(find . -name "*.bam" -exec basename {} \;)
do
	echo doing bamCoverage on $f
	bamCoverage -p 16 -b $f -o bw/coverage_$(cut -d'.' -f1 <<< $f).bw 
	echo coverage_$(cut -d'.' -f1 <<< $f).bw file created
done


## indexing not done
for f in $(find . -name "*.bam" -exec basename {} \;)
do
	echo "Indexing:"$f
	samtools index @ 16 $f
	echo $f".bai index file created"
	echo doing bamCoverage on $f
	bamCoverage -p 16 -b $f -o Coverage_$(cut -d'.' -f1 <<< $f).bw
	echo Coverage_$(cut -d'.' -f1 <<< $f).bw file created
done
