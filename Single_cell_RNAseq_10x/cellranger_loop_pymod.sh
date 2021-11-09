#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cellranger_loop
#BSUB -q bsc_debug
#BSUB -W 1:00
#BSUB -eo /gpfs/projects/cek26/SANTAMARIA_09_TMP/bsub_reports/cellranger_loop.err
#BSUB -oo /gpfs/projects/cek26/SANTAMARIA_09_TMP/bsub_reports/cellranger_loop.out
#BSUB -M 1800
#BSUB -n 1


# Load modules in order to run cellranger
module load singularity 
module load PYTHON/3.9.0

# Set and switch to your working directory
wd=/gpfs/projects/cek26/SANTAMARIA_09_TMP
cd $wd

# Set reference directory
ref=/gpfs/projects/cek26/genome_mus/cellranger/refdata-gex-mm10-2020-A

# Scripts directory
scr=/gpfs/projects/cek26/software/scripts

for f in $(find fastq_files -name  "*.gz")
do
echo $f
gzip -t $f && echo ok || echo corrupted
done

#Rename fastq files to cellranger input format
#Remember to setup the correct paths for the arguments
#To return help about the usage run: fastq_cellranger_format_converter.py -h

for f in $(find fastq_files -name  "*.gz")
do
echo $f
gzip -t $f && echo ok || echo corrupted
done
python $scr/fastq_cellranger_format_converter.py -t $wd/SANTAMARIA_09.txt -fd $wd/fastq_files -r False -cwd $wd
#Check if the files are renamed propperly and all expected are present
sort $wd/fastq_comp1.txt > $wd/sorted_fastq_comp1.txt
sort $wd/fastq_comp2.txt > $wd/sorted_fastq_comp2.txt
cmp $wd/sorted_fastq_comp1.txt $wd/sorted_fastq_comp2.txt || echo "Fastq files missing!"



# Create folders results and store FASTQfiles in a new folder [/fastq_files] inside the working directory, use the correct name nomenclature!
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# Where Read Type is one of:
# I1: Sample index read (optional)
# R1: Read 1
# R2: Read 2
mkdir cellranger_counts
mkdir to_bsub

# Loop to generate cellranger counts for each sample
# Firstly, generate a list of the samples you have
for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo $(cut -d'_' -f1 <<< $f) >> samples.txt
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

echo module load singularity
echo cd $wd/cellranger_counts
echo singularity exec /apps/CELLRANGER/SRC/images/cellranger-6.0.1.sif cellranger count --id=$f \
                 --transcriptome=$ref \
                 --fastqs=$wd/fastq_files \
                 --expect-cells=10000 \
                 --sample=$f
} > to_bsub/cellranger_count_$f.sh
sed -i -e 's/\r$//' to_bsub/cellranger_count_$f.sh
bsub < to_bsub/cellranger_count_$f.sh
done


# Next steps will be on hold until previous pipelines are over
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
singularity exec /apps/CELLRANGER/SRC/images/cellranger-6.0.1.sif cellranger aggr --id=$aggr_id \
                --csv=$wd/$aggr_id.csv \
                --normalize=mapped


# Remove used scripts
rm -r ../to_bsub
rm ../experiment.txt