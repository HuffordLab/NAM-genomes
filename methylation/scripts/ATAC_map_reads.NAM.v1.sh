#!/bin/bash
# Mapping the NAM ATAC-seq quality-processed reads to the NAM genomes

mapping_list=genome_mapping_combinations.txt
ATAC_home_dir=ATAC_output_parent_directory
atac_fastq_homedir=directory_containg_ATAC_fastq.gz_files
genome_home_dir=genome_parent_directory

cd ${home_dir}
while read LINE; do

  INPUT_GENOME=$(echo ${LINE} | cut -d ' ' -f1)
  RECIPIENT_GENOME=$(echo ${LINE} | cut -d ' ' -f2)
  out_dir=${ATAC_home_dir}/${INPUT_GENOME}/aligned
  index=${genome_home_dir}/${RECIPIENT_GENOME}/bowtie2_index/offrate4/${RECIPIENT_GENOME}.seed1_off4

  mkdir -p ${out_dir}
  cd ${out_dir}

  for REP in ${atac_fastq_homedir}/${INPUT_GENOME}.rep*.R1.fastq.gz; do

    rep_basename=`basename --suffix=.R1.fastq.gz "$REP"`
    fastq1=${ATAC_home_dir}/${INPUT_GENOME}/${rep_basename}.${RECIPIENT_GENOME}.R1.fastp.fastq
    fastq2=${ATAC_home_dir}/${INPUT_GENOME}/${rep_basename}.${RECIPIENT_GENOME}.R2.fastp.fastq

    OUT=${out_dir}/${rep_basename}.${RECIPIENT_GENOME}.align_reads.sh
    echo "#PBS -S /bin/bash" > ${OUT}
    echo "#PBS -q highmem_q" >> ${OUT}
    echo "#PBS -N ${rep_basename}.${RECIPIENT_GENOME}.align_ATAC_reads" >> ${OUT}
    echo "#PBS -l nodes=1:ppn=32" >> ${OUT}
    echo "#PBS -l mem=100gb" >> ${OUT}
    echo "#PBS -l walltime=200:00:00" >> ${OUT}
    echo "" >> ${OUT}
    echo "ml Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1" >> ${OUT}
    echo "" >> ${OUT}
    echo "cd ${out_dir}" >> ${OUT}
    echo "" >> ${OUT}
    echo "bowtie2 --local --very-sensitive-local --seed 1 -q --no-mixed --no-discordant --maxins 1000 -p \${PBS_NUM_PPN} \\" >> ${OUT}
    echo "-x ${index} -1 ${fastq1} -2 ${fastq2} -S ${out_dir}/${rep_basename}.${RECIPIENT_GENOME}.sam" >> ${OUT}
    qsub ${OUT}

  done
done < ${mapping_list}
