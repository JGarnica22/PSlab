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

# Switch to your working directory
cd /gpfs/projects/cek26/10x_BDC_INS

# Create folders results and store FASTQfiles in a new folder
mkdir cellranger_counts
mkdir to_bsub

# Loop to generate cellranger counts for each sample
# Firstly, generate a list of the sample you have
for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo $(cut -d'_' -f1 <<< $f)_$(cut -d'_' -f2 <<< $f) >> samples.txt
done

# Then generate the loop with the unique samples

for f in $(sort -u samples.txt)
do
{
echo \#!/bin/bash
echo \#BSUB cwd /gpfs/projects/cek26
echo \#BSUB -J cellranger_counts_$f
echo \#BSUB -q bsc_ls
echo \#BSUB -W 42:00
echo \#BSUB -eo /gpfs/projects/cek26/10x_BDC_INS/bsub_reports/cellranger_counts_$f.err
echo \#BSUB -oo /gpfs/projects/cek26/10x_BDC_INS/bsub_reports/cellranger_counts_$f.out
echo \#BSUB -M 3000
echo \#BSUB -n 256
echo \#BSUB "span[ptile=16]"
echo \#BSUB -x

echo module load CELLRANGER/3.1.0
echo cd /gpfs/projects/cek26/10x_BDC_INS/cellranger_counts
echo cellranger count --id=$f \
                 --transcriptome=/gpfs/projects/cek26/10x_BDC_INS/refdata-gex-mm10-2020-A \
                 --fastqs=/gpfs/projects/cek26/10x_BDC_INS/fastq_files \
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
echo $i,/gpfs/projects/cek26/10x_BDC_INS/cellranger_counts/$i/outs/molecule_info.h5 >> $aggr_id.csv
done


# Run cellranger aggr when counts finish
mkdir cellranger_aggr
cd cellranger_aggr
cellranger aggr --id=$aggr_id \
                --csv=/gpfs/projects/cek26/10x_BDC_INS/$aggr_id.csv \
                --normalize=mapped




# Remove used scripts
rm -r ../to_bsub
rm ../experiment.txt
