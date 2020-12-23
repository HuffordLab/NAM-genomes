#!/bin/bash
for fq in /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq/raw_data/RNASeq_plate_1-6/*fq.gz; do
nam=$(basename ${fq} | cut -f 1 -d "_")
ofq=$(basename ${fq} | cut -f 2- -d "_")
mkdir -p ${nam};
ln -s $fq $nam/${ofq};
done
