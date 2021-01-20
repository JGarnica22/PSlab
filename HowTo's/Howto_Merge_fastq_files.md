# How to merge fastq files :cat:

Sometimes you may have different fastq or fq files from the same sample or condition but resequenced or from different lanes or flowcells. As a consequence, you may want to merge them into a unique file to analyze it.

The best way is using `cat` UNIX tool (_concatenate_). There are different ways to do it: listing the files, by a pattern or through a loop:

### Listing:
````
cat file1.fq file2.fq file3.fq file4.fq file5.fq > merged.fq
````
### Pattern:
````
cat file*.fq > merged.fq
````
### Loop:
````
for file in NA24694_GCCAAT_L001_R${i}_*fastq.gz; do
        cat "$file" >> EA00694_GCCAAT_L001_R${i}.fastq.gz
    done
````


NOTE: `cat` also works for compressed files (`*.gz`)
