#!/bin/bash
# Prepares the Reference Genome for mapping as well as for using it with GATK pipeline
# You need to supply the referece genome as REF below or as:
# ./GATK_00_PrepareRef.sh your_genome.fasta
module load picard
module load samtools
module load bwa
module load bedtools2
module load bioawk
REF="$1"
#index genome for (a) picard, (b) samtools and (c) bwa
PICARD_HOME=$(dirname $(which picard))
# dict
java -Xmx100G -jar $PICARD_HOME/picard.jar CreateSequenceDictionary \
  REFERENCE=${REF} \
  OUTPUT=${REF%.*}.dict
# sam index
samtools faidx ${REF}
# bwa index
bwa index ${REF}
# windows (100kb)
bioawk -c fastx '{print $name"\t"length($seq)}' ${REF} > ${REF%.*}_length.txt
bedtools makewindows -w 100000 -g ${REF%.*}_length.txt > ${REF%.*}_100kb_coords.bed
java -Xmx100G -jar $PICARD_HOME/picard.jar BedToIntervalList \
  INPUT= ${REF%.*}_100kb_coords.bed \
  SEQUENCE_DICTIONARY= ${REF%.*}.dict \
  OUTPUT=${REF%.*}_100kb_gatk_intervals.list
