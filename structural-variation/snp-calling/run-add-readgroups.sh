#!/bin/bash
module load picard
PICARD_HOME="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina"
REF="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1.fasta"
ulimit -c unlimited

fastq=$(ls *.gz |head -n 1)
fline=$(zcat $fastq |head -n 1 |sed 's/^@//g')
file=$1
RGID=$(echo $fline |cut -f 3,4 -d ":")
RGSM=$(basename $(pwd) |cut -f 2- -d "_")
RGLB="${RGSM}-L001"
RGPU=001

echo -e "$RGID\t$RGSM\t$RGLB\t$RGPU"


java -Djava.io.tmpdir=$TMPDIR -Xmx50G -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
      I=${file} \
      O=${file%.*}_new.bam \
      RGID=$RGID \
      RGLB=$RGLB \
      RGPL=ILLUMINA \
      RGPU=$RGPU \
      RGSM=$RGSM

module load samtools
samtools index ${file%.*}_new.bam
