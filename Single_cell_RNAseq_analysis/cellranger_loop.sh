#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cellranger_loop
#BSUB -q sequential
#BSUB -W 72:00
#BSUB -eo /gpfs/projects/cek26/10x_BDC_INS/bsub_reports/cellranger_loop.err
#BSUB -oo /gpfs/projects/cek26/10x_BDC_INS/bsub_reports/cellranger_loop.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x

# Load modules in order to run cellranger
module load CELLRANGER/3.1.0

# Set and switch to your working directory
wd=/gpfs/projects/cek26/10x_BDC_INS
cd $wd

# Create folders results and store FASTQfiles in a new folder [/fastq_files] inside the working directory, use the correct name nomenclature!
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# Where Read Type is one of:
# I1: Sample index read (optional)
# R1: Read 1
# R2: Read 2
mkdir cellranger_counts
mkdir to_bsub

# Loop to generate cellranger counts for each sample
# Firstly, generate a list of the sample you have
for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo $(cut -d'_' -f1 <<< $f)_$(cut -d'_' -f2 <<< $f) >> samples.txt
done

# Then use a loop with the unique samples to generate respective jobs to be sent to bsub
# if expected number of cells differ, the jobs need to be changed manually
for f in $(sort -u samples.txt)
do
{
echo \#!/bin/bash
echo \#BSUB cwd $wd
echo \#BSUB -J cellranger_counts_$f
echo \#BSUB -q bsc_ls
echo \#BSUB -W 42:00
echo \#BSUB -eo $wd/bsub_reports/cellranger_counts_$f.err
echo \#BSUB -oo $wd/bsub_reports/cellranger_counts_$f.out
echo \#BSUB -M 3000
echo \#BSUB -n 256
echo \#BSUB "span[ptile=16]"
echo \#BSUB -x

echo module load CELLRANGER/3.1.0
echo cd $wd/cellranger_counts
echo cellranger count --id=$f \
                 --transcriptome=$wd/refdata-gex-mm10-2020-A \
                 --fastqs=$wd/fastq_files \
                 --expect-cells=10000 \
                 --sample=$f
} > to_bsub/cellranger_count_$f.sh
sed -i -e 's/\r$//' to_bsub/cellranger_count_$f.sh
bsub < to_bsub/cellranger_count_$f.sh
done


# Do not start next steps until previous pipelines are over
while [ $(sort -u samples.txt | wc -l) != $(find ./cellranger_counts -type f -name "*_info.h5" -print | wc -l) ]
do
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
echo $i,$wd/cellranger_counts/$i/outs/molecule_info.h5 >> $aggr_id.csv
done


# Run cellranger aggr when counts finish
mkdir cellranger_aggr
cd cellranger_aggr
cellranger aggr --id=$aggr_id \
                --csv=$wd/$aggr_id.csv \
                --normalize=mapped




# Remove used scripts
rm -r ../to_bsub
rm ../experiment.txt
