
## ATTENTION: PIPELINE NOT TESTED YET!!!

#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cr_arc_loop
#BSUB -q sequential
#BSUB -W 48:00
#BSUB -eo /gpfs/projects/cek26/project/bsub_reports/cr_arc_loop.err
#BSUB -oo /gpfs/projects/cek26/project/bsub_reports/cr_arc_loop.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x


###################################################################################################################################

## Edit paths of your project, if needed also change the path of .out and .err files ##
## it is assumed that you will create a folder (i.e project) to store the whole project and subfolders for your fastq files.

##################################################################################################################################

wd=/gpfs/projects/cek26/project #working directory

fastq_ATAC=/gpfs/projects/cek26/project/fastq_files/ATAC_fastq # path directory of fastq files with ATAC data
fastq_RNA=/gpfs/projects/cek26/project/fastq_files/GEX_fastq # path directory of fastq files with gene expression data
# Create folders results and store FASTQfiles in a new folder [i.e. /fastq_files] inside the working directory
# Use the correct name nomenclature!!
# 	SampleName_S1_L001_R1_001.fastq.gz
# 	Where Read Type is one of:
# 	I1: Sample index read (optional)
# 	R1: Read 1
# 	R2: Read 2

reference=/gpfs/projects/cek26/project/Ref # path to cellranger-arc-compatible reference file

id=PS # Tittle of the project

#################################################################################################################################

# Running pipeline, DO NOT EDIT!! Unless you want to change specifics

##################################################################################################################################
# Change to working directory
cd $wd

# Create needed directories
mkdir cellranger_arc to_bsub bsub_reports


# Loop to generate cellranger-arc count for each sample
# Firstly, generate a list of the sample you have

## TO DO: LIST OF SAMPLES PREFIX (ID)!!!!! samplex.txt or vector samples

# Then use a loop with the unique samples to generate respective jobs to be sent to bsub
# if expected number of cells differ, the jobs need to be changed manually

for f in $samples
do
	# Generate libraries csv
	{
		echo fastqs,sample,library_type
		echo $fastq_RNA,$f,Gene Expression
		echo $fastq_ATAC,$f,Chromatin Accessibility
	} > cellranger_arc/$f_libraries.csv

	{
		echo \#!/bin/bash
		echo \#BSUB cwd $wd
		echo \#BSUB -J cellranger_arc_$f
		echo \#BSUB -q bsc_ls
		echo \#BSUB -W 24:00
		echo \#BSUB -eo $wd/bsub_reports/cellranger_arc_$f.err
		echo \#BSUB -oo $wd/bsub_reports/cellranger_arc_$f.out
		echo \#BSUB -M 3000
		echo \#BSUB -n 256
		echo \#BSUB "span[ptile=16]"
		echo \#BSUB -x

		echo module load CELLRANGER-ARC/1.0 # check proper version!
		echo cd $wd/cellranger_arc
		echo cellranger-arc count --id=$f --reference=$reference --libraries=$f_libraries.csv
	} > to_bsub/cellranger_arc_$f.sh
	sed -i -e 's/\r$//' to_bsub/cellranger_arc_$f.sh
	bsub < to_bsub/cellranger_arc_$f.sh
done


# Do not start next steps until previous pipelines are over
while [ $(sort -u $samples | wc -l) != $(find ./cellranger_arc -type f -name "web_summary.html" -print | wc -l) ]
do
	echo cellranger-arc count still running	
	sleep 300;
done

echo PIPELINE COMPLETED, proceed to Seurat analysis!

