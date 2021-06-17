samtools merge -@ 6 iNKT_multiome_gex.bam treated/outs/gex_possorted_bam.bam control/outs/gex_possorted_bam.bam

samtools sort -@ 6 iNKT_multiome_gex.bam -o iNKT_multiome_gex_sorted.bam

samtools index -@ 6 iNKT_multiome_gex_sorted.bam

sinto fragments -p 7 -b iNKT_multiome_gex_sorted.bam -f gex_fragments.bed

sort -k1,1 -k2,2n gex_fragments.bed > gex_fragments.sort.bed
bgzip -@ 8 gex_fragments.sort.bed
tabix -p bed gex_fragments.sort.bed.gz


for i in $(find -name "Cluster*.bam" -exec basename {} \;)
do
echo sorting $i
samtools sort -@ 5 $i -o sorted_$i
echo indexing $i
samtools index -@ 5 sorted_$i
echo doing bigwigs files for $i
bamCoverage -p 5 -v -b sorted_$i -o cover_sorted_$i.bw
done




samtools sort -@ 5 iNKT_multiome_atac.bam -o iNKT_multiome_atac_sorted.bam

samtools index -@ 5 iNKT_multiome_atac_sorted.bam

sinto filterbarcodes -p 5 -b iNKT_multiome_gex_sorted.bam -c ../cond.txt


for i in $(find -name "Cluster*.bam" -exec basename {} \;)
do
sinto filterbarcodes -p 5 -b $i -c cond.txt
done


for i in $(find -name "merged_atac*.bam" -exec basename {} \;)
do
echo sorting $i
samtools sort -@ 5 $i -o sorted_$i
echo indexing $i
samtools index -@ 5 sorted_$i
echo doing bigwigs files for $i
bamCoverage -p 5 -b sorted_$i -o cover_sorted_$i.bw
done


for i in $(ls)
do
bsub < $i
done

sinto filterbarcodes -p 5 -b gex_possorted_bam.bam -c ../../cells_ctl.txt

cd ../..
for i in $(find ./treated -name "Cluster*.bam" -exec basename {} \;)
do
samtools merge -@ 6 merged_atac_$i ./treated/$i ./control/$i
done


for i in $(find -name "sorted_mer*.bam" -exec basename {} \;)
do
echo doing bigwigs files for $i
bamCoverage -p 4 -b $i -o cover_normBPM_$i.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads
done


bamCoverage -p 3 -b $i -o cover_normBPM_$i.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads

cd ..
bamCoverage -p 3 -b sorted_merged_atac_Cluster5.bam -o cover_normBPM_sorted_merged_atac_Cluster5.bam --binSize 20 --normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads
