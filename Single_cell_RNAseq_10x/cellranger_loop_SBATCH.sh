#!/bin/bash
#SBATCH -D /gpfs/projects/cek26/experiments
#SBATCH --job-name="Cellranger_inkt"
#SBATCH -q bsc_ls
#SBATCH -t 40:00:00
#SBATCH -e /gpfs/projects/cek26/experiments/SANTAMARIA_38/bsub_reports/cellranger_loop.err
#SBATCH -o /gpfs/projects/cek26/experiments/SANTAMARIA_38/bsub_reports/cellranger_loop.out
#SBATCH -n 16
#SBATCH --mail-type=all
#SBATCH --mail-user=garnica@clinic.cat


# Set and switch to your working directory
wd=/gpfs/projects/cek26/experiments/SANTAMARIA_38
cd $wd
fq=$wd/fastq_files
# --constraint=highmem

# Set reference directory
ref=/gpfs/projects/cek26/data/genomes/genome_mus/cellranger/refdata-gex-mm10-2020-A


# Create folders results and store FASTQfiles in a new folder [/fastq_files] inside the working directory, use the correct name nomenclature!
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# Where Read Type is one of:
# I1: Sample index read (optional)
# R1: Read 1
# R2: Read 2
cellranger_counts=cellranger_counts
cellranger_aggr=cellranger_aggr
mkdir $cellranger_counts to_bsub


# Loop to generate cellranger counts for each sample
# Firstly, generate a list of the samples you have
for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo "${f%%_S*}" >> samples.txt
done
f=$(sort -u samples.txt)

module load cellranger/6.1.2
cd $wd/cellranger_counts3
cellranger count --id=CD4mitet \
                 --transcriptome=$ref \
                 --fastqs=$fq \
                 --expect-cells=10000 \
                 --sample=$f \
                 --no-bam

# Then use a loop with the unique samples to generate respective jobs to be sent to bsub
# if expected number of cells differ, the jobs need to be changed manually
for f in $(sort -u samples.txt)
do
{
echo \#!/bin/bash
echo \#SBATCH -D $wd
echo \#SBATCH --job-name="CR_inkt_$f"
echo \#SBATCH -q bsc_ls
echo \#SBATCH -t 24:00:00
echo \#SBATCH -e $wd/bsub_reports/cellranger_loop_$f.err
echo \#SBATCH -e $wd/bsub_reports/cellranger_loop_$f.out
echo \#SBATCH -n 32
echo \#SBATCH --mail-type=all
echo \#SBATCH --mail-user=garnica@clinic.cat

echo module load cellranger/6.1.2
echo cd $wd/$cellranger_counts
echo cellranger count --id=$f \
                 --transcriptome=$ref \
                 --fastqs=$fq \
                 --expect-cells=10000 \
                 --sample=$f \
                 --no-bam
} > to_bsub/cellranger_count_$f.sh
sed -i -e 's/\r$//' to_bsub/cellranger_count_$f.sh
sbatch to_bsub/cellranger_count_$f.sh
done

# Next steps will be on hold until previous pipelines are over
while [ $(sort -u samples.txt | wc -l) != $(find ./cellranger_counts -type f -name "*_info.h5" -print | wc -l) ]
do
now=$(date)
echo Cellranger count still running... $now
sleep 300;
done

# Generate aggregation CSV 
for i in $(sort -u samples.txt)
do
echo $(cut -d'_' -f1 <<< $i) >> experiment.txt
done
aggr_id=()
for u in $(sort -u experiment.txt)
do
aggr_id=$aggr_id'_'$u
done
aggr_id=Aggr$aggr_id
echo library_id,molecule_h5 > $aggr_id.csv
for i in $(sort -u samples.txt)
do
echo $i,$wd/$cellranger_counts/$i/outs/molecule_info.h5 >> $aggr_id.csv
done


# Run cellranger aggr when counts finish
mkdir $cellranger_aggr
cd $cellranger_aggr
cellranger aggr --id=$aggr_id \
                --csv=$wd/$aggr_id.csv \
                --normalize=none




