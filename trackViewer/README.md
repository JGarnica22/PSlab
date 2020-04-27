# TrackViewer


## Prepare your input files:


#Install:

deeptools (bamCoverage - transform .bam to .bw)
samtools - index bam files (creates .bai)
bedtools - intersect genomic regions


Create bigwig from bam:

#1. Create index file (.bai) with samtools
````
samtools -b -o
````

* Remember to activate conda - you will see (base) in your terminal when conda is inititated
* If you need help on how to use samtools:
samtools --help


#2. Transform .bam to .bw with bamCoverage (deeptools)

* If you need help on how to use/what can be done with deeptools:
deeptools -h
or
deeptools --help
* Help on how to use bamCoverage:
bamCoverage --help


Use trackViewer: script.R
