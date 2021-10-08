#!/bin/bash

function check_data {
R1=$1;
R2=$(echo $R1 | sed ‘s/R1.fastq/R2.fastq/g’)

echo “Check if length of sequence equals quality score length”
bioawk -c fastx ‘{if (length($seq)!=length($qual)) print “Offending: ” NR}’ $R1
bioawk -c fastx ‘{if (length($seq)!=length($qual)) print “Offending: ” NR}’ $R2

echo “run FastQ Validator”
fastQValidator –file $R1
fastQValidator –file $R2
echo “Finished validation.”
echo “Check md5 sums:”;
echo `md5sum $R1`
echo `md5sum $R2`
echo “——–”
}

#validate Fastq file
echo “File: ” $1
check_data $1
