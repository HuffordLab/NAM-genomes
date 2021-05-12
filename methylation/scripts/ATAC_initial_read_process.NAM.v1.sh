#!/bin/bash

# The intial ATAC-seq read processing
# Workflow:
# (1) FASTP to trim and filter reads

mapping_list=genome_mapping_combinations
home_dir=ATAC_output_parent_directory
atac_fastq_homedir=directory_containg_ATAC_fastq.gz_files

while read LINE; do

  INPUT_GENOME=$(echo ${LINE} | cut -d ' ' -f1)
  RECIPIENT_GENOME=$(echo ${LINE} | cut -d ' ' -f2)
  fastp_cleaned_out=${home_dir}/${INPUT_GENOME}
  fastp_failed_out=${home_dir}/${INPUT_GENOME}
  cleaned_reports_out=${home_dir}/${INPUT_GENOME}
  out_dir=${home_dir}/${INPUT_GENOME}

  mkdir -p ${home_dir}/${INPUT_GENOME}
  cd ${home_dir}/${INPUT_GENOME}

  for REP in ${atac_fastq_homedir}/${INPUT_GENOME}.rep*.R1.fastq.gz; do
    rep_basename=`basename --suffix=.R1.fastq.gz "$REP"`
    in1=${atac_fastq_homedir}/${rep_basename}.R1.fastq.gz
    in2=${atac_fastq_homedir}/${rep_basename}.R2.fastq.gz

    OUT=${home_dir}/${INPUT_GENOME}/${rep_basename}.${RECIPIENT_GENOME}.process_reads.sh
    echo "#PBS -S /bin/bash" > ${OUT}
    echo "#PBS -q highmem_q" >> ${OUT}
    echo "#PBS -N ${rep_basename}.${RECIPIENT_GENOME}.ATAC_seq_read_process" >> ${OUT}
    echo "#PBS -l nodes=1:ppn=9" >> ${OUT}
    echo "#PBS -l mem=6gb" >> ${OUT}
    echo "#PBS -l walltime=8:00:00" >> ${OUT}
    echo "# Mapping ${rep_basename} to ${RECIPIENT_GENOME}" >> ${OUT}
    echo "" >> ${OUT}
    echo "threads=\$((\$PBS_NUM_PPN-5))" >> ${OUT}
    echo "" >> ${OUT}
    echo "cd ${out_dir}" >> ${OUT}
    echo "fastp \\" >> ${OUT}
    echo "--in1 ${in1} \\" >> ${OUT}
    echo "--in2 ${in2} \\" >> ${OUT}
    echo "--out1 ${fastp_cleaned_out}/${rep_basename}.${RECIPIENT_GENOME}.R1.fastp.fastq \\" >> ${OUT}
    echo "--out2 ${fastp_cleaned_out}/${rep_basename}.${RECIPIENT_GENOME}.R2.fastp.fastq \\" >> ${OUT}
    echo "--unpaired1 ${fastp_failed_out}/${rep_basename}.${RECIPIENT_GENOME}.R1.unpaired.fastp.fastq \\" >> ${OUT}
    echo "--unpaired2 ${fastp_failed_out}/${rep_basename}.${RECIPIENT_GENOME}.R2.unpaired.fastp.fastq \\" >> ${OUT}
    echo "--failed_out ${fastp_failed_out}/${rep_basename}.${RECIPIENT_GENOME}.failed.fastp.fastq \\" >> ${OUT}
    echo "--reads_to_process 0 \\" >> ${OUT}
    echo "--dont_overwrite \\" >> ${OUT}
    echo "--detect_adapter_for_pe \\" >> ${OUT}
    echo "--correction \\" >> ${OUT}
    echo "--length_required 35 \\" >> ${OUT}
    echo "--json ${cleaned_reports_out}/${rep_basename}.${RECIPIENT_GENOME}.fastp.json \\" >> ${OUT}
    echo "--html ${cleaned_reports_out}/${rep_basename}.${RECIPIENT_GENOME}.fastp.html \\" >> ${OUT}
    echo "--report_title \"${rep_basename} mapped to ${RECIPIENT_GENOME}\" \\" >> ${OUT}
    echo "--thread \${threads} \\" >> ${OUT}
    echo "2>${cleaned_reports_out}/${rep_basename}.${RECIPIENT_GENOME}.fastp.summary_report" >> ${OUT}
    qsub ${OUT}

  done
done < ${mapping_list}
