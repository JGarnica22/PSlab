#!/bin/bash
#SBATCH -D /gpfs/projects/cek26
#SBATCH --job-name="hello_test"
#SBATCH -q debug
#SBATCH -t 2:00:00
#SBATCH -e /gpfs/projects/cek26/experiments/test.err
#SBATCH -o /gpfs/projects/cek26/experiments/test.out
#SBATCH -n 16
#SBATCH --mail-type=all
#SBATCH --mail-user=jmorobor@gmail.com


echo Trying Slurm
