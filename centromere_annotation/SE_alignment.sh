#!/bin/bash

chip_SE="$1"
genome="$2"
prefix=$(basename $chip_SE |cut -f1 -d ".")

module load BWA/0.7.17-GCC-8.3.0
cd /scratch/jl03308/NAM_pancentromere/analysis/peak_call/${prefix}
index=${genome%.%*}.ann
if [ -f "${index}" ]; then
   echo "skip genome indexing."
else
  bwa index $genome
fi
bwa mem -t 12 ${genome} $chip_SE > ${prefix}.ChIP.sam

ml SAMtools/1.9-GCC-8.3.0
samtools view -@ 12 -b -o ${prefix}.ChIP.bam ${prefix}.ChIP.sam
samtools sort -o ${prefix}.ChIP.sorted.bam -T ${prefix}.ChIP -@ 12 ${prefix}.ChIP.bam

samtools view -@ 12 -bhq 20 ${prefix}.ChIP.sorted.bam -o ${prefix}.ChIP.q20.bam
samtools sort -o ${prefix}.ChIP.q20.sorted.bam -T ${prefix}.ChIP -@ 12 ${prefix}.ChIP.q20.bam
samtools flagstat ${prefix}.ChIP.q20.sorted.bam > ${prefix}.ChIP.flagstat
samtools index -@ 12 ${prefix}.ChIP.q20.sorted.bam
