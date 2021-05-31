

#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cr_arc_loop_24
#BSUB -q sequential
#BSUB -W 72:00
#BSUB -eo /gpfs/projects/cek26/SANTAMARIA_23_24/bsub_reports/cr_arc_loop_24.err
#BSUB -oo /gpfs/projects/cek26/SANTAMARIA_23_24/bsub_reports/cr_arc_loop_24.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x


###################################################################################################################################

## Edit paths of your project, if needed also change the path of .out and .err files ##
## it is assumed that you will create a folder (i.e project) to store the whole project and subfolders for your fastq files.

##################################################################################################################################

wd=/gpfs/projects/cek26/SANTAMARIA_23_24 #working directory

fastq_ATAC=/gpfs/projects/cek26/SANTAMARIA_23_24/fastq_files/ATAC_fastq # path directory of fastq files with ATAC data
# ATAC data should have index fastq files (I2)!! Files from CNAG are UMI2.
fastq_RNA=/gpfs/projects/cek26/SANTAMARIA_23_24/fastq_files/GEX_fastq # path directory of fastq files with gene expression data
# Create folders results and store FASTQfiles in a new folder [i.e. /fastq_files] inside the working directory
# Use the correct name nomenclature!!
# 	SampleName_S1_L001_R1_001.fastq.gz
# 	Where Read Type is one of:
# 	I1: Sample index read (optional)
# 	R1: Read 1
# 	R2: Read 2

reference=/gpfs/projects/cek26/genome_mus/cellranger_arc/refdata-cellranger-arc-mm10-2020-A-2.0.0 # path to cellranger-arc-compatible reference file

id=TR1 # Tittle of the project

#################################################################################################################################

# Running pipeline, DO NOT EDIT!! Unless you want to change specifics

##################################################################################################################################
# Change to working directory
cd $wd

# Create needed directories
mkdir cellrangerarc_counts


# Loop to generate cellranger-arc count for each sample
# Firstly, generate a list of the sample you have

for f in $(find ./fastq_files -name "*.fastq.gz" -exec basename {} \;)
do
echo $(cut -d'_' -f1 <<< $f) >> samples.txt
done

# Then use a loop with the unique samples to generate respective jobs to be sent to bsub
# if expected number of cells differ, the jobs need to be changed manually

for f in $(sort -u samples.txt)
do
	# Generate libraries csv
	{
		echo fastqs,sample,library_type
		echo $fastq_RNA,$f"_GEX",Gene Expression
		echo $fastq_ATAC,$f"_ATAC",Chromatin Accessibility
	} > cellrangerarc_counts/$f"_libraries".csv

	{
		echo \#!/bin/bash
		echo \#BSUB cwd $wd
		echo \#BSUB -J cellranger_arc_$f
		echo \#BSUB -q bsc_ls
		echo \#BSUB -W 48:00
		echo \#BSUB -eo $wd/bsub_reports/cellranger_arc_$f.err
		echo \#BSUB -oo $wd/bsub_reports/cellranger_arc_$f.out
		echo \#BSUB -M 3000
		echo \#BSUB -n 256
		echo \#BSUB "span[ptile=16]"
		echo \#BSUB -x

		echo module purge
    echo module load singularity # check proper version!
		echo cd $wd/cellrangerarc_counts
		echo singularity exec /gpfs/apps/MN3/CELLRANGER-ARC/SRC/images/cellranger-arc-2.0.0.sif cellranger-arc count --id=$f --reference=$reference --libraries=$f"_libraries".csv
	} > to_bsub/cellranger_arc_$f.sh
	sed -i -e 's/\r$//' to_bsub/cellranger_arc_$f.sh
	bsub < to_bsub/cellranger_arc_$f.sh
done


# Do not start next steps until previous pipelines are over
while [ $(sort -u samples.txt | wc -l) != $(find ./cellrangerarc_counts -type f -name "cloupe.cloupe" -print | wc -l) ]
do
	echo cellranger-arc count still running	
	sleep 500;
done



# Generate aggregation CSV 

echo library_id,atac_fragments,per_barcode_metrics,gex_molecule_info > aggr_$id.csv
for i in $(sort -u samples.txt)
do
echo $i,$wd/cellrangerarc_counts/$i/outs/atac_fragments.tsv.gz,$wd/cellrangerarc_counts/$i/outs/per_barcode_metrics.csv,$wd/cellrangerarc_counts/$i/outs/gex_molecule_info.h5 >> aggr_$id.csv
done


# Run cellranger aggr when counts finish
module purge
module load singularity
mkdir cellranger_arc_aggr
cd cellranger_arc_aggr
singularity exec /gpfs/apps/MN3/CELLRANGER-ARC/SRC/images/cellranger-arc-2.0.0.sif cellranger-arc aggr --id=$id \
                --csv=$wd/aggr_$id.csv \
                --normalize=none \
                --reference=$reference



  echo PIPELINE $id COMPLETED, proceed to Seurat analysis!




