#!/bin/bash
nam=$1
ls *.sortedByCoord.out.bam > input.fofn
module load samtools
samtools merge --threads 36 -O BAM -b input.fofn merged_${nam}.bam
