#!/bin/bash
module load star
genome="$1"
index="$(basename ${genome%.*})-star"
base=$(basename $(pwd))
read1=$2
read2=$(echo $read1 |sed 's/_R1.fq.gz/_R2.fq.gz/g')
out=$(basename ${read1%%.*} | sed 's/_R1.fq.gz/_/g')

STAR \
--genomeDir $index \
--runThreadN 36 \
--sjdbFileChrStartEnd SJ.all \
--runMode alignReads \
--readFilesIn $read1 $read2 \
--readFilesCommand zcat \
--outSAMattributes All \
--outSAMmapqUnique 10 \
--outFilterMismatchNmax 0 \
--outFileNamePrefix ${out}_round-2 \
--outBAMsortingThreadN 4 \
--limitBAMsortRAM 5594394835 \
--outSAMtype BAM SortedByCoordinate \
--outWigType bedGraph read1_5p


#chmod go+rwx -R
#module load fastqc
#fastqc $read1
#fastqc $read2
