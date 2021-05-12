#!/bin/bash
# submit_process_NAM_genomes file

genome_list=/scratch/war88452/NAM_founders/genomes/NAM_founder_genome_list.txt
while read GENOME; do
  qsub -v genome=${GENOME} process_NAM_genomes.sh
done < <(cut -f1 ${genome_list} | sort -u)
