#!/bin/bash
module load star
genome="$1"
index="$(basename ${genome%.*})-star"
base=$(basename $(pwd))
read1="$2"
read2=$(echo $read1 |sed 's/_R1.fq.gz/_R2.fq.gz/g')
out=$(basename ${read1%%.*} | sed 's/_R1.fq.gz/_/g')


STAR \
--genomeDir $index \
--runThreadN 36 \
--runMode alignReads \
--readFilesIn $read1 $read2 \
--readFilesCommand zcat \
--outFileNamePrefix ${out}_ \
--outSAMtype None
