#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J HMMRATAC
#BSUB -q bsc_ls
#BSUB -W 16:00
#BSUB -eo /gpfs/projects/cek26/ATACseq/bsub_reports/HMMRATAC.err
#BSUB -oo /gpfs/projects/cek26/ATACseq/bsub_reports/HMMRATAC.out
#BSUB -M 3000
#BSUB -n 128
#BSUB -R "span[ptile=16]"
#BSUB -x

# Load necessary modules
module purge
module load java/1.8.0u66 samtools 

# Set your working directory
wd=/gpfs/projects/cek26/ATACseq
cd $wd

# For peak calling with HMMRATAC you only can use pair-end data
# Before starting, you need to index your bam files into .bai and generate your genome.info
# Easiest way to do it, do this command on one of your bam files batch:
# samtools view -H *.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info
# Run HMMRATAC
mkdir peak_calling/HMMRATAC
for f in $(find ./alignment -name "trim_bowtie*.bam" -exec basename {} \;)
do
java -jar ../software/HMMRATAC_V1.2.10_exe.jar -b alignment/$f -i alignment/$f.bai \
-g alignment/genome.info -o peak_calling/HMMRATAC/HMM_
done
