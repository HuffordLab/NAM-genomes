#!/bin/bash

# Extract_all_Cs_NAM_genomes.v4.sh
# This version extracts the cytosine sequence contexts from the NAM founder genomes
# It works on the 26 genomes in parallel

home_dir=genome_parent_directory
genome_list=NAM_founder_genome_list.txt

cd ${home_dir}

while read genome; do

  mkdir -p ${home_dir}/${genome}/split_by_chromosome
  mkdir -p ${home_dir}/${genome}/C_coordinates_all_contexts

  cd ${home_dir}/${genome}

  OUT=${home_dir}/${genome}/${genome}.extract_Cs.sh
  echo "#PBS -S /bin/bash" >> ${OUT}
  echo "#PBS -q highmem_q" >> ${OUT}
  echo "#PBS -N extract_Cs_${genome}" >> ${OUT}
  echo "#PBS -l nodes=1:ppn=8" >> ${OUT}
  echo "#PBS -l mem=100gb" >> ${OUT}
  echo "#PBS -l walltime=14:00:00" >> ${OUT}
  echo "" >> ${OUT}

  echo "ml pyfaidx/0.5.5.1-foss-2016b-Python-2.7.14" >> ${OUT}
  echo "ml EMBOSS/6.6.0-foss-2016b" >> ${OUT}

  # The chromosomes are split into separate files because that's the only way fuzznuc runs
  echo "cd ${home_dir}/${genome}/split_by_chromosome" >> ${OUT}
  echo "faidx --split-files ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.fa" >> ${OUT}

    # Run fuzznuc on the split fasta files by chromosome for CG, CHG, and CHH
  echo "out_dir=${home_dir}/${genome}/C_coordinates_all_contexts" >> ${OUT}
  echo "cd \${out_dir}" >> ${OUT}

  echo "for context in cgn chg chh; do" >> ${OUT}
  echo "  for chrom in ${home_dir}/${genome}/split_by_chromosome/*.fa ; do" >> ${OUT}
  echo "    basename=\`basename --suffix=.fa \"\$chrom\"\`" >> ${OUT}
  echo "    fuzznuc \${chrom} -pattern \${context} -complement -rformat gff -outfile \${out_dir}/\${basename}.\${context}.gff" >> ${OUT}
  echo "  done" >> ${OUT}

  # Now concatenate the chromosomes and sort by coordinates
  echo "  cat \${out_dir}/*.\${context}.gff | grep -v '^#' \\" >> ${OUT}
  echo "  | awk -v var=\"\$context\" '{if (\$7 == \"+\") {print \$1,\$4-1,\$4,var,\".\",\"+\"} else if (\$7 == \"-\") {print \$1,\$5-1,\$5,var,\".\",\"-\"}}' OFS=\"\t\" \\" >> ${OUT}
  echo "  > \${out_dir}/concatenated.\${context}.bed" >> ${OUT}

  echo "  sort -b -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} \${out_dir}/concatenated.\${context}.bed > \${out_dir}/${genome}.\${context}.sorted.bed" >> ${OUT}
  echo "  awk '{print \$1\":\"\$2\"..\"\$3,\$1,\$2,\$3,\$6}' OFS=\"\t\" \${out_dir}/${genome}.\${context}.sorted.bed > \${out_dir}/${genome}.\${context}.sorted.to_join_temp1.bed" >> ${OUT}
  echo "  sort -b -k1,1 --parallel=\${PBS_NUM_PPN} \${out_dir}/${genome}.\${context}.sorted.to_join_temp1.bed > \${out_dir}/${genome}.\${context}.sorted.to_join.bed" >> ${OUT}

  echo "done" >> ${OUT}
  qsub ${OUT}

done < <(cut -f1 ${genome_list} | sort -u)
