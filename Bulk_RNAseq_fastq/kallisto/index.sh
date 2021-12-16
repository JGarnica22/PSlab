#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --qos=debug
#SBATCH --time=02:00:00
#SBATCH --error=/gpfs/projects/cek26/data/genomes/genome_mus/kallisto/kallisto_index.err
#SBATCH --output=/gpfs/projects/cek26/data/genomes/genome_mus/kallisto/kallisto_index.out
#SBATCH -N 4
#SBATCH --cpus-per-task=16



# Load necessary modules
module load intel/2021.4.0 impi/2021.4.0 hdf5 szip kallisto

# Set your working directory
wd=/gpfs/projects/cek26/data/genomes/genome_mus/kallisto
cd $wd

kallisto index -i GRCm38.p6.idx ../GRCm38.p6.genome.fa
