#!/bin/bash
#SBATCH --job-name=kallisto_loop
#SBATCH --qos=debug
#SBATCH --time=00:20:00
#SBATCH --error=/gpfs/projects/cek26/experiments/kallisto/bsub_reports/kallisto_loop.err
#SBATCH --output=/gpfs/projects/cek26/experiments/kallisto/bsub_reports/kallisto_loop.out


# Set your working directory
wd=/gpfs/projects/cek26/experiments/kallisto
gdir=/gpfs/projects/cek26/data/genomes/genome_mus
to_sbatch=to_sbatch
cd $wd
mkdir alignment $to_sbatch

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
{
echo \#!/bin/bash
echo \#SBATCH --job-name=$(cut -d'_' -f2 <<< $f)
echo \#SBATCH --qos=debug
echo \#SBATCH --time=02:00:00
echo \#SBATCH --error=$wd/bsub_reports/$(cut -d'_' -f2 <<< $f).err
echo \#SBATCH --output=$wd/bsub_reports/$(cut -d'_' -f2 <<< $f).out
echo \#SBATCH -N 1
echo \#SBATCH -n 2
echo \#SBATCH -c 1

echo module load intel/2021.4.0 impi/2021.4.0 hdf5 szip kallisto
echo kallisto quant -i $gdir/kallisto/GRCm38.p6.idx -o alignment -b 50 --single -t 16 -l 200 -s 30 \<\(zcat $f\)

} > $to_sbatch/align_k_$(cut -d'_' -f2 <<< $f).sh
sed -i -e 's/\r$//' $to_sbatch/align_k_$(cut -d'_' -f2 <<< $f).sh
sbatch $to_sbatch/align_k_$(cut -d'_' -f2 <<< $f).sh
done




