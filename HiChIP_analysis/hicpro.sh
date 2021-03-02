#!/bin/bash
#BSUB -cwd /gpfs/projects/cek26
#BSUB -J hicpro
#BSUB -q sequential
#BSUB -W 48:00
#BSUB -eo /gpfs/projects/cek26/HiC-Pro/bsub_reports/hicpro.err
#BSUB -oo /gpfs/projects/cek26/HiC-Pro/bsub_reports/hicpro.out
#BSUB -M 3000
#BSUB -n 16
#BSUB -R "span[ptile=16]"

module load singularity

# Switch to your working directory
wd=/gpfs/projects/cek26/HiC-Pro
cd $wd

singularity exec /apps/HIC-PRO/SRC/images/HiC-Pro-3.0.0.sif /HiC-Pro_3.0.0/bin/HiC-Pro -i /gpfs/projects/cek26/HiC-Pro/data -o /gpfs/projects/cek26/HiC-Pro/out -c /gpfs/projects/cek26/HiC-Pro/config-hicpro.txt