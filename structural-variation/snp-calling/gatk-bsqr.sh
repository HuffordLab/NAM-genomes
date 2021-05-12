#!/bin/bash
module load gatk
module load r-geneplotter
REF=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1.fasta
GATK="java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar /work/LAS/mhufford-lab/arnstrm/programs/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar"
FRVCF="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/c-gatk/gatk-first-round-files/first_round_sorted-snps-only_filtered-pass-only.vcf"
IBAM="$1"
OBAM=$(basename ${IBAM%.*})

${GATK} BaseRecalibrator \
   --reference $REF \
   --input $IBAM \
   --known-sites ${FRVCF} \
   --output ${OBAM}_bef-R1.table

${GATK} ApplyBQSR \
   --reference $REF \
   --input $IBAM \
   --output ${OBAM}_recal.bam \
   --bqsr-recal-file ${OBAM}_bef-R1.table

${GATK} BaseRecalibrator \
   --reference $REF \
   --input ${OBAM}_recal.bam \
   --known-sites ${FRVCF} \
   --output ${OBAM}_aft-R1.table

${GATK} AnalyzeCovariates \
   -before ${OBAM}_bef-R1.table \
   -after ${OBAM}_aft-R1.table \
   -plots ${OBAM}-AnalyzeCovariates.pdf
