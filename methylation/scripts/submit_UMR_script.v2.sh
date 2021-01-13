#!/bin/bash

# Running qsubs for calling UMRs in the NAM methylomes

OUTDIR=UMR_output_directory
merged_list=NAM_methylome_list_merge.v2.txt

while read LINE; do

  PREFIX=$(echo ${LINE} | cut -d ' ' -f1)
  GENOME_NAME=$(echo ${LINE} | cut -d ' ' -f2)

  qsub -v prefix=${PREFIX},out_dir=${OUTDIR},GENOME=${GENOME_NAME} UMR_algorithm_all_contexts_v1.3.NAM_founders.sh

done < <(grep -v 'skip' ${merged_list} | cut -f1,4 | sort -u)
