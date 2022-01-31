#!/bin/bash
#SBATCH -D /gpfs/projects/cek26
#SBATCH --job-name="hello_test"
#SBATCH -q bsc_ls
#SBATCH -t 48:00:00
#SBATCH -e /gpfs/projects/cek26/experiments/test.err
#SBATCH -o /gpfs/projects/cek26/experiments/test.out
#SBATCH -n 64
#SBATCH --mail-type=all
#SBATCH --mail-user=jmorobor@gmail.com


echo Trying Slurm
