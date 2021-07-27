wd=/gpfs/projects/cek26/SANTAMARIA_15/tracer


cd $wd

#add file name to the fasta format to have all sequences labeled
ls > ../cells.txt

for i in $(sort -u ../cells.txt)
do
cd $i/unfiltered_TCR_seqs
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' $i"_TCRseqs.fa" > $i"_TCRseqs_ti.fa"
cd ../..
done


# Merge all fasta files

find . -type f -name '*ti.fa' -exec cat {} ';' > merged_all.fasta


# upload this file to http://www.imgt.org/HighV-QUEST/search.action
