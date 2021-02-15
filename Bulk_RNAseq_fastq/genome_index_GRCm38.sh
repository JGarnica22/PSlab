#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J STAR_index
#BSUB -q bsc_debug
#BSUB -W 02:00
#BSUB -eo /gpfs/projects/cek26/mm10/reports/STAR_index.err
#BSUB -oo /gpfs/projects/cek26/mm10/reports/STAR_index.out
#BSUB -M 1800

# Load modules in order to run FastQC and STAR
module purge
module load java/1.8.0u66 fastqc intel/2017.4 impi/2017.4 mkl/2017.4 gcc/5.3.0 gcc/4.9.1 openssl/1.1.1c python/3.7.4_pip STAR/2.7.5a

# Create a STAR index for the genome we are going to align to

# Download and then upload to cluster your genome of reference and respective GTF annotation file. Both files 
# should be downloaded from same source GENCODE or Ensembl to avoid nomenclature problems. To download from GENCODE use:

# Get both files (current release) from https://www.gencodegenes.org/mouse.html and decompress both files by:
# gzip -d <>

# Create index for your genome of reference and annotation file
cd mm10

STAR --runMode genomeGenerate --genomeDir /gpfs/projects/cek26/mm10/STAR_index --genomeFastaFiles GRCm38.p6.genome.fa --sjdbGTFfile gencode.vM25.annotation.gtf --sjdbOverhang 49 --runThreadN 12
