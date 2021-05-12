#!/bin/bash

# filename:
# convert_CGmap_to_allc_then_merge_then_convert_bed.v2.sh

# workflow:
# a separate script is submitted for each separate 'TISSUE'. There are 51 total 'TISSUES':
  # 1. convert CGmap files to allc files. The names of the CGmap files are in field 5 of the NAM_methylome_list_merge.v2.txt file
  # 2. merge the allc files
  # 3. convert the allc files to bed files
  # 4. process the bed files to get coverage, bigwigs, etc, for viewing in genome browser

# software:
# methylpy 1.3.2
# BEDTools 2.28.0

out_dir=output_directory
CGmap_home=directory_containing_CGmap_files
merge_list=NAM_methylome_list_merge.v2.txt
home_dir=genome_parent_directory

cd ${out_dir}

while read LIB_LINE; do
  TISSUE=$(echo ${LIB_LINE} | cut -d ' ' -f1)
  GENOME=$(echo ${LIB_LINE} | cut -d ' ' -f2)

  OUT=${out_dir}/${TISSUE}.sh
  echo "#PBS -S /bin/bash" > ${OUT}
  echo "#PBS -q highmem_q" >> ${OUT}
  echo "#PBS -N ${TISSUE}.convertCGmap_allc" >> ${OUT}
  echo "#PBS -l nodes=1:ppn=4" >> ${OUT}
  echo "#PBS -l mem=32gb" >> ${OUT}
  echo "#PBS -l walltime=24:00:00" >> ${OUT}
  echo "" >> ${OUT}
  echo "# This script was produced by ${out_dir}/$0" >> ${OUT}
  echo "# This script was produced on `date`" >> ${OUT}
  echo "" >> ${OUT}
  echo "module load methylpy/1.3.2-foss-2016b-Python-2.7.14" >> ${OUT}
  echo "ml BEDTools/2.28.0-foss-2018b" >> ${OUT}
  echo "cd ${out_dir}" >> ${OUT}

# convert CGmap to allc

  while read REP; do
    basename=`basename --suffix=.CGmap "${REP}"`
    echo "awk '{if ( \$2 == \"C\" && \$4 != \"--\") {print \$1,\$3,\"+\",\$4,\$7,\$8,\"1\"} \\" >> ${OUT}
    echo "else if  ( \$2 == \"G\" && \$4 != \"--\") {print \$1,\$3,\"-\",\$4,\$7,\$8,\"1\"} }' OFS=\"\t\" ${CGmap_home}/${REP} > ${out_dir}/${basename}.uncorrected.allc" >> ${OUT}
# correct scaffold names to match the reference genome
    echo "sed 's/scaf_alt_/scaf-alt_/' ${out_dir}/${basename}.uncorrected.allc > ${out_dir}/${basename}.allc" >> ${OUT}
    echo "rm ${out_dir}/${basename}.uncorrected.allc" >> ${OUT}
  done < <(awk -v TISSUE="$TISSUE" '{if ($1==TISSUE) {print $5}}' OFS="\t" ${merge_list} | cut -f1 | sort | uniq)

# merge the allc reps

  echo "methylpy merge-allc --allc-files \\" >> ${OUT}
  while read REP; do
    echo "${out_dir}/${REP} \\" >> ${OUT}
  done < <(awk -v TISSUE="$TISSUE" '{if ($1==TISSUE) {print $2}}' OFS="\t" ${merge_list} | cut -f1 | sort | uniq)
  echo "--output-file ${out_dir}/${TISSUE}.merged_reps.allc --num-procs 12 --compress-output False" >> ${OUT}

# final step: process the coverages into bed files, split by sequence contexts, and convert to bigwigs
# workflow:
# (1) convert the allc file to a bed file
# (2) split the bed file by sequence contexts. This is modified for the NAM contexts, which contain Hs
# (3) Create the files containing coverage for all Cs in the genome, split into the three contexts, and also CHGCGN
# (4) Create coverage files for the three contexts

  echo "" >> ${OUT}
  allc=${out_dir}/${TISSUE}.merged_reps.allc

  # convert allc to allc.bed and create a new file
  echo "awk '{print \$1,\$2-1,\$2,\$3,\$4,\$5,\$6}' OFS=\"\\t\" ${allc} > ${out_dir}/${TISSUE}.allc.bed.cnn" >> ${OUT}
  # extract the CG context Cs from the new bed file and create a new sorted file
  echo "awk '{if (\$5==\"CG\") {print \$0}}' OFS=\"\\t\" ${out_dir}/${TISSUE}.allc.bed.cnn > ${out_dir}/${TISSUE}.cgn.temp1" >> ${OUT}
  echo "sort -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.cgn.temp1 > ${out_dir}/${TISSUE}.allc.bed.cgn" >> ${OUT}
  # extract the CHG context Cs from the new bed file and create a new sorted file
  echo "awk '{if (\$5==\"CHG\") {print \$0}}' OFS=\"\\t\" ${out_dir}/${TISSUE}.allc.bed.cnn > ${out_dir}/${TISSUE}.chg.temp1" >> ${OUT}
  echo "sort -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.chg.temp1 > ${out_dir}/${TISSUE}.allc.bed.chg" >> ${OUT}
  # extract the CHH context Cs from the new bed file and create a new sorted file
  echo "awk '{if (\$5==\"CHH\") {print \$0}}' OFS=\"\\t\" ${out_dir}/${TISSUE}.allc.bed.cnn > ${out_dir}/${TISSUE}.chh.temp1" >> ${OUT}
  echo "sort -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.chh.temp1 > ${out_dir}/${TISSUE}.allc.bed.chh" >> ${OUT}
  echo "rm ${out_dir}/${TISSUE}.cgn.temp1 ${out_dir}/${TISSUE}.chg.temp1 ${out_dir}/${TISSUE}.chh.temp1" >> ${OUT}
  # combine the allc files for CGN and CHG to create a CHGCGN allc bed file
  echo "cat ${out_dir}/${TISSUE}.allc.bed.chg ${out_dir}/${TISSUE}.allc.bed.cgn \\" >> ${OUT}
  echo "| sort -b -k1,1 -k2,2n -k3,3n > ${out_dir}/${TISSUE}.allc.bed.chgcgn" >> ${OUT}

  # Create the files of total Cs in the genome with coverage and methylation information
  # Doing this for each of the three sequence contexts
  echo "for context in cgn chg chh; do" >> ${OUT}
  # Prepare to join the allc bed file to the total C bed file (for each sequence context)
  echo "awk '{print \$1\":\"\$2\"..\"\$3,\$1,\$2,\$3,\$6,\$7}' OFS=\"\\t\" ${out_dir}/${TISSUE}.allc.bed.\${context} > ${out_dir}/${TISSUE}.allc.bed.\${context}.to_join.temp1" >> ${OUT}
  echo "sort -b -k1,1 --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.allc.bed.\${context}.to_join.temp1 > ${out_dir}/${TISSUE}.allc.bed.\${context}.to_join" >> ${OUT}
  # Join
  echo "join -a 1 -e 0 -o 1.2,1.3,1.4,1.5,2.5,2.6 -1 1 -2 1 \\" >> ${OUT}
  echo "${home_dir}/${GENOME}/C_coordinates_all_contexts/${GENOME}.\${context}.sorted.to_join.bed \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.allc.bed.\${context}.to_join -t \$'\\t' \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage.bed.unsorted" >> ${OUT}
  echo "awk '{ if ((\$6 > 0) && (\$5/\$6) < 0.2 ) {print \$0,(\$5/\$6),\"less\"} \\" >> ${OUT}
  echo "else if ((\$6 > 0) && (\$5/\$6) >= 0.2 ) {print \$0,(\$5/\$6),\"greater\"} \\" >> ${OUT}
  echo "else if (\$6 == 0) {print \$0,\"2\",\"greater\" }}' OFS=\"\\\t\" ${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage.bed.unsorted \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage_ratios.bed.unsorted" >> ${OUT}
  echo "sort -b -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage_ratios.bed.unsorted \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.total_\${context}s_coverage_ratios.bed" >> ${OUT}
  # Remove temp files
  echo "rm \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.allc.bed.\${context}.to_join.temp1 \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.allc.bed.\${context}.to_join \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage.bed.unsorted \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.allc.bed.total_\${context}s_coverage_ratios.bed.unsorted" >> ${OUT}
  echo "done" >> ${OUT}
  # combine the CGN and CHG total files
  echo "cat ${out_dir}/${TISSUE}.total_chgs_coverage_ratios.bed ${out_dir}/${TISSUE}.total_cgns_coverage_ratios.bed \\" >> ${OUT}
  echo "| sort -b -k1,1 -k2,2n -k3,3n > ${out_dir}/${TISSUE}.total_chgcgns_coverage_ratios.bed" >> ${OUT}

  # create the bedgraph files for the three sequence contexts, % methylation and coverage
  echo "for context in cgn chg chh; do" >> ${OUT}
  echo "awk '{print \$1,\$2,\$3,\$7}' OFS=\"\t\" ${out_dir}/${TISSUE}.allc.bed.\${context} \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.coverage.unsorted.bg.\${context}" >> ${OUT}
  echo "sort -b -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.coverage.unsorted.bg.\${context} \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.coverage.sorted.bg.\${context}" >> ${OUT}
  echo "awk '{print \$1,\$2,\$3,(\$6/\$7)}' OFS=\"\t\" ${out_dir}/${TISSUE}.allc.bed.\${context} \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.methylation.unsorted.bg.\${context}" >> ${OUT}
  echo "sort -b -k1,1 -k2,2n -k3,3n --parallel=\${PBS_NUM_PPN} ${out_dir}/${TISSUE}.methylation.unsorted.bg.\${context} \\" >> ${OUT}
  echo "> ${out_dir}/${TISSUE}.methylation.sorted.bg.\${context}" >> ${OUT}

  # make bigwigs for coverage and methylation files
  echo "/home/war88452/local_software/kent_tools/bedGraphToBigWig \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.coverage.sorted.bg.\${context} \\" >> ${OUT}
  echo "${home_dir}/${GENOME}/processed_genome/${GENOME}.reduced_header.bedtools_genome_file.tsv \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.coverage.sorted.\${context}.bw" >> ${OUT}
  echo "/home/war88452/local_software/kent_tools/bedGraphToBigWig \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.methylation.sorted.bg.\${context} \\" >> ${OUT}
  echo "${home_dir}/${GENOME}/processed_genome/${GENOME}.reduced_header.bedtools_genome_file.tsv \\" >> ${OUT}
  echo "${out_dir}/${TISSUE}.methylation.sorted.\${context}.bw" >> ${OUT}
  echo "rm ${out_dir}/${TISSUE}.coverage.unsorted.bg.\${context}" >> ${OUT}
  echo "rm ${out_dir}/${TISSUE}.coverage.sorted.bg.\${context}" >> ${OUT}
  echo "rm ${out_dir}/${TISSUE}.methylation.unsorted.bg.\${context}" >> ${OUT}
  echo "rm ${out_dir}/${TISSUE}.methylation.sorted.bg.\${context}" >> ${OUT}
  echo "done" >> ${OUT}

  qsub ${OUT}

done < <(awk '{if ($3=="merge") {print $0}}' OFS="\t" ${merge_list} | cut -f1,4 | sort | uniq)
