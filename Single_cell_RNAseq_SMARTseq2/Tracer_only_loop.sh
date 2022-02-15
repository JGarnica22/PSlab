#!/bin/bash
#SBATCH -D /gpfs/projects/cek26/experiments
#SBATCH --job-name="TRACER"
#SBATCH -q bsc_ls
#SBATCH -t 48:00:00
#SBATCH -e /gpfs/projects/cek26/experiments/SANTAMARIA_34/bsub_reports/tracer.err
#SBATCH -o /gpfs/projects/cek26/experiments/SANTAMARIA_34/bsub_reports/tracer.out
#SBATCH -n 4
#SBATCH --mail-type=all
#SBATCH --mail-user=garnica@clinic.cat


wd=/gpfs/projects/cek26/experiments/SANTAMARIA_34
cd $wd

tracer=tracer
mkdir $tracer to_bsub

# Split fastq files in groups to parallelize the jobs, here we will be splitting the files in 3:
part=$(expr $(ls fastq_files/*_1.fastq.gz | wc -l) / 3)

for i in \$Group_1 \$Group_2 \$Group_3
do

# Loop for TraCer
{
echo \#!/bin/bash
echo \#SBATCH -D /gpfs/projects/cek26/experiments
echo \#SBATCH --job-name="TRACER"
echo \#SBATCH -q bsc_ls
echo \#SBATCH -t 48:00:00
echo \#SBATCH -e /gpfs/projects/cek26/experiments/SANTAMARIA_34/bsub_reports/$i.err
echo \#SBATCH -o /gpfs/projects/cek26/experiments/SANTAMARIA_34/bsub_reports/$i.out
echo \#SBATCH -n 64
echo \#SBATCH --mail-type=all
echo \#SBATCH --mail-user=garnica@clinic.cat

echo module load singularity

echo cd $wd
echo all=\$\(find ./fastq_files -name "*_1.fastq.gz" -exec basename {} "\;" \)
echo Group_1\=\$\(echo \"\$all\" \| head \-n $part\)
echo Group_2\=\$\(echo \"\$all\" \| tail \-n $part\)
echo Group_3\=\$\(echo \"\$all\" \| tail \-n $(expr $part \* 2 + 2) \| head \-n $(expr $part + 2)\)

# Loop for TraCer
echo for f in $i 
echo do
echo if \[ \-d $tracer\/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) \]\; then
echo echo \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) "done"
echo else
echo echo doing tracer \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)
echo singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif tracer assemble -p 64 -s Hsap --loci A B --max_junc_len 50 \
fastq_files/\$f fastq_files/\$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\)_2.fastq.gz \$\(cut \-d \".\" \-f1 \<\<\< \$f \| sed \'s/..$//\'\) $tracer
echo fi
echo done
} > to_bsub/tracer_$i.sh
sed -i -e 's/\r$//' to_bsub/tracer_$i.sh
sbatch to_bsub/tracer_$i.sh
done
 

# Summarise TraCer
while [ $(ls $tracer| wc -l) != $(find ./fastq_files -type f -name "*_1.fastq.gz" -print | wc -l) ]
do
echo Tracer still running
sleep 300;
done


module load singularity
echo running Tracer summarise
singularity exec /apps/TRACER/SRC/images/tracer-0.6.0.sif tracer summarise -p 4 -s Mmus $tracer
