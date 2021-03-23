#!/bin/bash
gff=$1
bam=$2
ml purge
ml subread
featureCounts -T 36 -a $gff -o ${gff%.*}_.counts.txt -t mRNA -g Name -p --tmpDir $TMPDIR  -B -C --primary bam-files/${bam}_*.sortedByCoord.out.bam

