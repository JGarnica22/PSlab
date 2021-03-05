
## ATTENTION: PIPELINE NOT TESTED YET!!!

#!/bin/bash
#BSUB cwd /gpfs/projects/cek26
#BSUB -J cr_multi_loop
#BSUB -q sequential
#BSUB -W 48:00
#BSUB -eo /gpfs/projects/cek26/project/bsub_reports/cr_multi_loop.err
#BSUB -oo /gpfs/projects/cek26/project/bsub_reports/cr_multi_loop.out
#BSUB -M 1800
#BSUB -n 16
#BSUB -R "span[ptile=16]"
#BSUB -x


###################################################################################################################################

## Edit paths of your project, if needed also change the path of .out and .err files ##

##################################################################################################################################

wd=/gpfs/projects/cek26/project #working directory

fastq_AC=/gpfs/projects/cek26/project/fastq_files # path directory of fastq files with features or ADT data
fastq_RNA=/gpfs/projects/cek26/project/fastq_files # path directory of fastq files with gene expression data
fastq_VDJ=/gpfs/projects/cek26/project/fastq_files # path directory of fastq files with v(d)j data
# Create folders results and store FASTQfiles in a new folder [i.e. /fastq_files] inside the working directory
# Use the correct name nomenclature!!
# 	[Sample Name-type(FB/GEX/VDJ)]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# 	Where Read Type is one of:
# 	I1: Sample index read (optional)
# 	R1: Read 1
# 	R2: Read 2

reference_RNA=  # Path of folder containing 10x-compatible reference.
reference_AC= #/path/to/feature/reference, Feature reference CSV file, declaring Feature Barcode constructs and associated barcodes.
reference_VDJ= # Path of folder containing 10x-compatible VDJ reference. Required for Immune Profiling libraries.

id=PS # Tittle of the project
exp_cells=3000 #Expected number of recovered cells. Default: 3,000 cells. Switch code into --force-cells to force number of cells,
			   # bypassing the cell detection algorithm. Use this if the number of cells estimated by Cell Ranger is 
			   # not consistent with the barcode rank plot.
bam=true #skip generation of bam file, if false .bam and .bai files will be generated, if true not.

## IMPORTANT: cellranger aggr is set for samples coming from differents donors (see donor and origin at 10X_multi_VDJ_RNA_ADT.md)
## if otherwise aggr csv config file must be modified!!



#################################################################################################################################

# Running pipeline, DO NOT EDIT!! Unless you want to change specifics

##################################################################################################################################
# Change to working directory
cd $wd

# Create needed directories
mkdir cellranger_multi to_bsub bsub_reports


# Loop to generate cellranger multi for each sample
# Firstly, generate a list of the sample you have

## TO DO: LIST OF SAMPLES PREFIX (ID)!!!!! samplex.txt or vector samples

# Then use a loop with the unique samples to generate respective jobs to be sent to bsub
# if expected number of cells differ, the jobs need to be changed manually

for f in $samples
do
	# Generate multi config csv
	{
		echo \# $f multi config file,,,
		echo [gene-expression],,,
		echo reference,$reference_RNA,,
		echo expect-cells,$exp_cells,,
		echo no-bam,$bam,,
		echo ,,,
		echo [feature],,,
		echo reference,$reference_AC,,
		echo ,,,
		echo [vdj],,,
		echo reference,$reference_VDJ
		echo ,,,
		echo [libraries],,,
		echo fastq_id,fastqs,lanes,feature_types
		echo $f-FB,$fastq_AC,any,Antibody Capture
		echo $f-FB,$fastq_RNA,any,Gene Expression
		echo $f-FB,$fastq_VDJ,any,VDJ
	} > cellranger_multi/$f_multi_config.csv

	{
		echo \#!/bin/bash
		echo \#BSUB cwd $wd
		echo \#BSUB -J cellranger_multi_$f
		echo \#BSUB -q bsc_ls
		echo \#BSUB -W 24:00
		echo \#BSUB -eo $wd/bsub_reports/cellranger_multi_$f.err
		echo \#BSUB -oo $wd/bsub_reports/cellranger_multi_$f.out
		echo \#BSUB -M 3000
		echo \#BSUB -n 256
		echo \#BSUB "span[ptile=16]"
		echo \#BSUB -x

		echo module load CELLRANGER/6.0
		echo cd $wd/cellranger_multi
		echo cellranger multi --id=$f --csv=$f_multi_config.csv
	} > to_bsub/cellranger_multi_$f.sh
	sed -i -e 's/\r$//' to_bsub/cellranger_multi_$f.sh
	bsub < to_bsub/cellranger_multi_$f.sh
done


# Do not start next steps until previous pipelines are over
while [ $(sort -u $samples | wc -l) != $(find ./cellranger_multi -type f -name "web_summary.html" -print | wc -l) ]
do
	echo cellranger multi still running	
	sleep 300;
done

echo generating aggregation CSV

# Generate aggregation CSV
echo library_id,library_outs,donor,origin > aggr_id.csv 
for f in $samples
do
echo $f,$wd/cellranger_multi/$f/outs/per_sample_outs/$f,$f,$f >> aggr_id.csv
done

echo running cellranger aggr
module load CELLRANGER/6.0 #check proper version!!
# Run cellranger aggr
mkdir cellranger_multi_aggr && cd cellranger_multi_aggr

cellranger aggr --id=$id \
                --csv=$wd/aggr_id.csv \
                --normalize=none


echo cellranger aggr completed
echo PIPELINE $id COMPLETED, proceed to Seurat analysis!
