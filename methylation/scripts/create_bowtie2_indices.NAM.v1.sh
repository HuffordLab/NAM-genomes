#!/bin/bash

# This creates scripts to make the Bowtie2 indices for the 25 NAM founders and B73 v5
# It works on the 26 genomes in parallel

home_dir=genome_parent_directory
genome_list=NAM_founder_genome_list.txt

cd ${home_dir}

while read genome; do

  mkdir -p ${home_dir}/${genome}/bowtie2_index/offrate4
  cd ${home_dir}/${genome}
  OUT=${home_dir}/${genome}/${genome}.bowtie2_index.sh
  echo "#PBS -S /bin/bash" >> ${OUT}
  echo "#PBS -q highmem_q" >> ${OUT}
  echo "#PBS -N ${genome}.bowtie2_index" >> ${OUT}
  echo "#PBS -l nodes=1:ppn=6" >> ${OUT}
  echo "#PBS -l mem=20gb" >> ${OUT}
  echo "#PBS -l walltime=12:00:00" >> ${OUT}
  echo "" >> ${OUT}
  echo "ml Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1" >> ${OUT}
  echo "cd ${home_dir}/${genome}/bowtie2_index/offrate4" >> ${OUT}
  echo "bowtie2-build --seed 1 --offrate 4 --threads 6 \\" >> ${OUT}
  echo "${home_dir}/${genome}/processed_genome/${genome}.reduced_header.fa \\" >> ${OUT}
  echo "${home_dir}/${genome}/bowtie2_index/offrate4/${genome}.seed1_off4" >> ${OUT}
  qsub ${OUT}

done < <(cut -f1 ${genome_list} | sort -u)
