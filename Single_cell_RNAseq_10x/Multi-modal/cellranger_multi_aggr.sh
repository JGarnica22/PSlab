#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cellranger_aggr_Parvus
#BSUB -q sequential
#BSUB -W 72:00
#BSUB -eo /gpfs/projects/cek26/Parvus_10X/bsub_reports/cellranger_aggr_Parvus.err
#BSUB -oo /gpfs/projects/cek26/Parvus_10X/bsub_reports/cellranger_aggr_Parvus.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x

# Load modules in order to run cellranger
module load CELLRANGER/5.0

# Set and switch to your working directory
wd=/gpfs/projects/cek26/Parvus_10X
cd $wd


# echo library_id,library_outs,donor,origin > aggr_id.csv
for i in $(find cellranger -maxdepth 1 -name  "*multi*.csv" -exec basename {} \;)
do
echo $(cut -d '_' -f1 <<< $i),$wd/cellranger/$(cut -d '_' -f1 <<< $i)/outs,$(cut -d '_' -f1 <<< $i),$(cut -d '-' -f1 <<< $i) >> aggr_id.csv
done



# Run cellranger aggr when counts finish
mkdir cellranger_aggr
cd cellranger_aggr
cellranger aggr --id=Parvus \
                --csv=$wd/aggr_id.csv \
                --normalize=mapped
