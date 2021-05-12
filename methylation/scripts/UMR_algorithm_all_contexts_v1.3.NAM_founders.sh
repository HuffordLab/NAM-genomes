#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N UMR_algorithm_all_contexts_v1.3.NAM_founders
#PBS -l nodes=1:ppn=4
#PBS -l mem=200gb
#PBS -l walltime=12:00:00

# This is a batch-submission version with a separate batch-submission script.
# The $GENOME, $out_dir, and $prefix variables are submitted externally
# The separate script is called submit_UMR_script.v2.sh
# This version is designed to be run on the 26 NAM founders, a total 51 'TISSUES'


ml BEDTools/2.28.0-foss-2018b

home_dir=genome_parent_directory
allc_bed_home=directory_containing_merged_allc_files_from_step_six
genome=${home_dir}/${GENOME}/processed_genome/${GENOME}.reduced_header.bedtools_genome_file.tsv
w100_s20=${home_dir}/${GENOME}/processed_genome/${GENOME}.w100_s20.bed
contigs_no_genes=${home_dir}/${GENOME}/processed_genome/${GENOME}.chromosomes_contigs_no_genes.bed
blacklist=${home_dir}/${GENOME}/processed_genome/${GENOME}.reduced_header.N_gaps.bed
total_genes=${home_dir}/${GENOME}/processed_genome/${GENOME}.genes.bed

mkdir -p ${out_dir}
cd ${out_dir}

# for the loop below, running through the whole algorithm with both chg and chgcgn inputs
for context in chg chgcgn; do

  allc_context=${allc_bed_home}/${prefix}.allc.bed.${context}
  total_Cs_context=${allc_bed_home}/${prefix}.total_${context}s_coverage_ratios.bed

  # Create a new allc.bed CHG file but for only CHGs with ratio < 0.2
  awk '{if ( ($6/$7) < 0.2 ) {print $0}}' OFS="\t" ${allc_context} > ${out_dir}/${prefix}.allc.${context}_less0.2.bed

  # Create the initial merged UMRs
  bedtools intersect -wa -wb -a ${w100_s20} -b ${allc_context} > ${out_dir}/${prefix}.${context}.temp1
  sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} ${out_dir}/${prefix}.${context}.temp1 > ${out_dir}/${prefix}.${context}.temp2
  bedtools groupby -i ${out_dir}/${prefix}.${context}.temp2 -g 1,2,3 -c 9,10,1 -o sum,sum,count > ${out_dir}/${prefix}.${context}.temp3
  awk '{if ( $5 >= 20 && $6 >= 4 && ( ($4/$5) < 0.2 ) ) {print $0}}' OFS="\t" ${out_dir}/${prefix}.${context}.temp3 > ${out_dir}/${prefix}.${context}.temp4
  bedtools intersect -v -wa -a ${out_dir}/${prefix}.${context}.temp4 -b ${blacklist} \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} \
  | bedtools merge -d 30 -i stdin > ${out_dir}/${prefix}.${context}.temp5

  # Trim the UMR edges
  bedtools intersect -wa -wb -a ${out_dir}/${prefix}.${context}.temp5 -b ${out_dir}/${prefix}.allc.${context}_less0.2.bed > ${out_dir}/${prefix}.${context}.temp6
  sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} ${out_dir}/${prefix}.${context}.temp6 > ${out_dir}/${prefix}.${context}.temp7
  bedtools groupby -i ${out_dir}/${prefix}.${context}.temp7 -g 1,2,3 -c 6,6 -o min,max > ${out_dir}/${prefix}.${context}.temp8
  awk '{print $1,$4-1,$5}' OFS="\t" ${out_dir}/${prefix}.${context}.temp8 > ${out_dir}/${prefix}.${context}.temp9

  # Remove UMRs with a total mC/cov ratio >= 20, coverage < 20, and counts < 4
  bedtools intersect -wa -wb -a ${out_dir}/${prefix}.${context}.temp9 -b ${allc_context} > ${out_dir}/${prefix}.${context}.temp10
  sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} ${out_dir}/${prefix}.${context}.temp10 > ${out_dir}/${prefix}.${context}.temp11
  bedtools groupby -i ${out_dir}/${prefix}.${context}.temp11 -g 1,2,3 -c 9,10,1 -o sum,sum,count > ${out_dir}/${prefix}.${context}.temp12
  awk '{ if ( ( ($4/$5) < 0.2 ) && ( $5 >= 20 ) && ( $6 >= 4 ) ) {print $0,(100*($4/$5))}}' OFS="\t" ${out_dir}/${prefix}.${context}.temp12 > ${out_dir}/${prefix}.w100_s20.UMRs.no_fill-ins.${context}.bed

  #############################################################################################################################
  #############################################################################################################################

  # Create the UMR complement and split by less or greater than 4kb in length
  bedtools complement -i ${out_dir}/${prefix}.w100_s20.UMRs.no_fill-ins.${context}.bed -g ${genome} > ${out_dir}/${prefix}.complement_total.${context}.bed
  awk '{ if (($3-$2) <= 4000) {print $1,$2+1,$3}}' OFS="\t" ${out_dir}/${prefix}.complement_total.${context}.bed > ${out_dir}/${prefix}.complement.${context}.bed
  awk '{ if (($3-$2) > 4000) {print $1,$2+1,$3}}' OFS="\t" ${out_dir}/${prefix}.complement_total.${context}.bed > ${out_dir}/${prefix}.complement_greater_4k.${context}.bed

  #############################################################################################################################
  #############################################################################################################################

  # Intersect the complements with the total CHG coverage file
  bedtools intersect -wa -wb -loj -a ${out_dir}/${prefix}.complement.${context}.bed -b ${total_Cs_context} \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.intersect_total_${context}.bed

  # Complements lacking any CHGs are set aside for merging ("merge_comps1.bed").
  awk '{ if ($11 == ".") {print $1,$2,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.intersect_total_${context}.bed > ${out_dir}/${prefix}.merge_comps1.${context}.bed

  # Collapse the complement-CHG intersects with groupby
  bedtools groupby -i ${out_dir}/${prefix}.complement.intersect_total_${context}.bed -g 1,2,3 -c 9,9,10,10 -o min,max,min,max \
  > ${out_dir}/${prefix}.complement.intersect_total_${context}.grouped_cov_ratios.bed

  # Complements for which all CHGs have coverage and a per nucleotide mC/cov ratio < 0.2 are set aside for merging ("merge_comps2.bed").
  awk '{ if (($4 > 0) && ($7 < 0.2)) {print $1,$2,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.intersect_total_${context}.grouped_cov_ratios.bed \
  > ${out_dir}/${prefix}.merge_comps2.${context}.bed

  # Prepare the grouped complements for intersect parsing
  awk '{ if (($6 < 0.2) && ( ! (($4 > 0) && ($7 < 0.2))) ) {print $1,$2,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.intersect_total_${context}.grouped_cov_ratios.bed \
  > ${out_dir}/${prefix}.use_for_parsing.${context}.bed

  # Intersect the filtered complements with the total CHG coverage file
  bedtools intersect -wa -wb -loj -a ${out_dir}/${prefix}.use_for_parsing.${context}.bed -b ${total_Cs_context} \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.bed

  # Separate the intersected complements into "less" and "greater" CHGs
  awk '{ if ($11 == "less") {print $0}}' OFS="\t" ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.bed > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.bed
  awk '{ if ($11 == "greater") {print $0}}' OFS="\t" ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.bed > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.bed

  # Collapse with groupby using the interval for both the "less" and "greater" file
  bedtools groupby -i ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.bed -g 1,2,3 -c 6,6 -o min,max > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed
  bedtools groupby -i ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.bed -g 1,2,3 -c 6,6 -o min,max > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed

  # prep the grouped less and greater files for joining together
  awk '{print $1":"$2".."$3,$0}' OFS="\t" ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed | sort -b -k1,1 --parallel=${PBS_NUM_PPN} \
  > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed_to_join
  awk '{print $1":"$2".."$3,$0}' OFS="\t" ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed | sort -b -k1,1 --parallel=${PBS_NUM_PPN} \
  > ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed_to_join

  # Join them together with the field order: less min coord, less max coord, greater min coord, greater max coord
  join -a 1 -e 0 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -1 1 -2 1 \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed_to_join \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed_to_join \
  -t $'\t' > ${out_dir}/${prefix}.${context}.final_parse

  # Separate the suitable left and right edges of the complements
  awk '{ if ($4 < $6) {print $1,$2,$6}}' OFS="\t" ${out_dir}/${prefix}.${context}.final_parse > ${out_dir}/${prefix}.${context}.final_parse.1.left
  awk '{ if ($5 > $7) {print $1,$7,$3}}' OFS="\t" ${out_dir}/${prefix}.${context}.final_parse > ${out_dir}/${prefix}.${context}.final_parse.1.right

  # Intersect the edges of the complements with the CHGs with per nucleotide ratios < 0.2
  bedtools intersect -wa -wb -a ${out_dir}/${prefix}.${context}.final_parse.1.left -b ${out_dir}/${prefix}.allc.${context}_less0.2.bed \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.${context}.final_parse.1.left.intersect
  bedtools intersect -wa -wb -a ${out_dir}/${prefix}.${context}.final_parse.1.right -b ${out_dir}/${prefix}.allc.${context}_less0.2.bed \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.${context}.final_parse.1.right.intersect

  # The final filter to create "merge_comps3.bed" and "merge_comps4.bed"
  bedtools groupby -i ${out_dir}/${prefix}.${context}.final_parse.1.left.intersect -g 1,2,3 -c 6 -o max \
  | awk '{print $1,$2,$4}' OFS="\t" > ${out_dir}/${prefix}.merge_comps3.${context}.bed

  bedtools groupby -i ${out_dir}/${prefix}.${context}.final_parse.1.right.intersect -g 1,2,3 -c 6 -o min \
  | awk '{print $1,$4,$3}' OFS="\t" > ${out_dir}/${prefix}.merge_comps4.${context}.bed

  #############################################################################################################################
  #############################################################################################################################

  # For complements > 4kb, take the left and right inner flanks of 2 kb (called 'L_IF_2kb' and 'R_IF_2kb')
  awk '{print $1,$2,($2+2000),$3}' OFS="\t" ${out_dir}/${prefix}.complement_greater_4k.${context}.bed \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.L_IF_2kb.${context}.bed

  awk '{print $1,($3-2000),$3,$2}' OFS="\t" ${out_dir}/${prefix}.complement_greater_4k.${context}.bed \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.R_IF_2kb.${context}.bed

  # Intersect the L_IF_2kb and R_IF_2kb complements with the total CHG coverage file
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools intersect -wa -wb -loj -a ${out_dir}/${prefix}.complement.${IF_2kb}.${context}.bed -b ${total_Cs_context} \
    | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.bed
  done

  # Collapse the complement-CHG intersects with groupby.
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools groupby -i ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.bed -g 1,2,3,4 -c 10,10,11,11 -o min,max,min,max \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.grouped_cov_ratios.bed
  done

  #############################################################################################################################
  #############################################################################################################################

  # The complement inner flanks for which all CHGs have coverage and a per nucleotide mC/cov ratio < 0.2 are expanded back to original complements.
  awk '{ if (($5 > 0) && ($8 < 0.2)) {print $1,$2,$4}}' OFS="\t" ${out_dir}/${prefix}.complement.L_IF_2kb.intersect_total_${context}.grouped_cov_ratios.bed \
  > ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}.bed

  awk '{ if (($5 > 0) && ($8 < 0.2)) {print $1,$4,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.R_IF_2kb.intersect_total_${context}.grouped_cov_ratios.bed \
  > ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}.bed

  # Intersect the expanded complements with the total CHG coverage file
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools intersect -wa -wb -loj -a ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}.bed -b ${total_Cs_context} \
    | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}_intersect.bed
  done

  # Filter the intersected expanded complements for just "greater" CHGs
  awk '{ if ($11 == "greater") {print $0}}' OFS="\t" ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.bed > ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.bed
  awk '{ if ($11 == "greater") {print $0}}' OFS="\t" ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.bed > ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.bed

  # Collapse the intersected expanded complements and set the inner boundaries with a greater CHG coordinate
  bedtools groupby -i ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.bed -g 1,2,3 -c 6 -o min > ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.bed
  bedtools groupby -i ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.bed -g 1,2,3 -c 6 -o max > ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.bed

  awk '{print $1,$2,$4}' OFS="\t" ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.bed \
  > ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed

  awk '{print $1,$4,$3}' OFS="\t" ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.bed \
  > ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed

  # Intersect these reduced inner flanks with the CHG coverage file containing only CHGs < 0.2 ratio
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools intersect -wa -wb -a ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}_intersect.greater.coord.reduced.bed \
    -b ${out_dir}/${prefix}.allc.${context}_less0.2.bed \
    | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}_intersect.greater.coord.reduced.bed.intersect_less_0.2.bed
  done

  # The final grouping and filter to prep the expanded inner flanks for merging to the UMRs
  bedtools groupby -i ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed.intersect_less_0.2.bed \
  -g 1,2,3 -c 6 -o max | awk '{print $1,$2,$4}' OFS="\t" \
  > ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.to_merge.${context}.bed

  bedtools groupby -i ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed.intersect_less_0.2.bed \
  -g 1,2,3 -c 6 -o min | awk '{print $1,$4,$3}' OFS="\t" \
  > ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.to_merge.${context}.bed

  #############################################################################################################################
  #############################################################################################################################

  # Now to deal with L_IF_2kb and R_IF_2kb that contain both CHGs < 0.2 and CHGs >= 0.2 and lacking coverage

  # Filter out the inner flanks that were already processed. Keep only inner flanks that have >= 1 lesser CHG
  # and >= 1 CHG that is greater or lacking coverage.
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    awk '{ if (($7 < 0.2) && ($8 >= 0.2)) {print $1,$2,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.grouped_cov_ratios.bed \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.${context}.bed
  done

  # Intersect the mixed L_IF_2kb and R_IF_2kb with the total CHG coverage file
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools intersect -wa -wb -loj -a ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.${context}.bed -b ${total_Cs_context} \
    | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.bed
  done

  # Separate the intersected inner flanks into "less" and "greater" CHGs
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    awk '{ if ($11 == "less") {print $0}}' OFS="\t" ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.bed \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.bed
    awk '{ if ($11 == "greater") {print $0}}' OFS="\t" ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.bed \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.bed
  done

  # Collapse with groupby using the interval for both the "less" and "greater" file
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools groupby -i ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.bed -g 1,2,3 -c 6,6 -o min,max \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed
    bedtools groupby -i ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.bed -g 1,2,3 -c 6,6 -o min,max \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed
  done

  # prep the grouped less and greater files for joining together
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    awk '{print $1":"$2".."$3,$0}' OFS="\t" ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed | sort -b -k1,1 --parallel=${PBS_NUM_PPN} \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed.to_join
    awk '{print $1":"$2".."$3,$0}' OFS="\t" ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed | sort -b -k1,1 --parallel=${PBS_NUM_PPN} \
    > ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed.to_join
  done

  # Join them together with the field order: 4:  min(Less coord), 5: max(Less coord), 6: min(Greater coord), 7: max(Greater coord)
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    join -a 1 -e 0 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -1 1 -2 1 \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed.to_join \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed.to_join \
    -t $'\t' > ${out_dir}/${prefix}.complement.${IF_2kb}.to_parse.${context}.bed
  done

  # Filter for inner flanks where the lesser CHG is outside of the greater CHG, then set the new inside coordinate
  awk '{ if ($4 < $6) {print $1,$2,$6}}' OFS="\t" ${out_dir}/${prefix}.complement.L_IF_2kb.to_parse.${context}.bed > ${out_dir}/${prefix}.complement.L_IF_2kb.filtered.${context}.bed
  awk '{ if ($5 > $7) {print $1,$7,$3}}' OFS="\t" ${out_dir}/${prefix}.complement.R_IF_2kb.to_parse.${context}.bed > ${out_dir}/${prefix}.complement.R_IF_2kb.filtered.${context}.bed

  # Intersect the filtered inner flanks with the CHGs with per nucleotide ratios < 0.2
  for IF_2kb in L_IF_2kb R_IF_2kb; do
    bedtools intersect -wa -wb -a ${out_dir}/${prefix}.complement.${IF_2kb}.filtered.${context}.bed -b ${out_dir}/${prefix}.allc.${context}_less0.2.bed \
    | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} > ${out_dir}/${prefix}.complement.${IF_2kb}.filtered.intersect.${context}.bed
  done

  # Set the final inner coordinates for the filtered inner flanks. The less CHG sets the coordinate
  bedtools groupby -i ${out_dir}/${prefix}.complement.L_IF_2kb.filtered.intersect.${context}.bed -g 1,2,3 -c 6 -o max \
  | awk '{print $1,$2,$4}' OFS="\t" > ${out_dir}/${prefix}.complement.L_IF_2kb.to_merge.${context}.bed

  bedtools groupby -i ${out_dir}/${prefix}.complement.R_IF_2kb.filtered.intersect.${context}.bed -g 1,2,3 -c 6 -o min \
  | awk '{print $1,$4,$3}' OFS="\t" > ${out_dir}/${prefix}.complement.R_IF_2kb.to_merge.${context}.bed

  #############################################################################################################################
  #############################################################################################################################

  # Prior to the final merge, must filter all merge files to remove all blacklist overlaps. Catting together then filtering.
  # The final merge: nine bed files are being merged together
  cat \
  ${out_dir}/${prefix}.w100_s20.UMRs.no_fill-ins.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.to_merge.${context}.bed \
  ${out_dir}/${prefix}.merge_comps1.${context}.bed \
  ${out_dir}/${prefix}.merge_comps2.${context}.bed \
  ${out_dir}/${prefix}.merge_comps3.${context}.bed \
  ${out_dir}/${prefix}.merge_comps4.${context}.bed \
  | awk '{print $1,$2,$3}' OFS="\t" \
  | bedtools intersect -v -wa -a stdin -b ${blacklist} \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} \
  | bedtools merge -d 5 -i stdin > ${out_dir}/${prefix}.UMRs_final_merge.${context}.bed

  #############################################################################################################################
  #############################################################################################################################
  # Remove the temp files

  for i in {1..12}; do
    rm ${out_dir}/${prefix}.${context}.temp${i}
  done

  rm \
  ${out_dir}/${prefix}.w100_s20.UMRs.no_fill-ins.${context}.bed \
  ${out_dir}/${prefix}.complement_total.${context}.bed \
  ${out_dir}/${prefix}.complement.${context}.bed \
  ${out_dir}/${prefix}.complement_greater_4k.${context}.bed \
  ${out_dir}/${prefix}.complement.intersect_total_${context}.bed \
  ${out_dir}/${prefix}.merge_comps1.${context}.bed \
  ${out_dir}/${prefix}.complement.intersect_total_${context}.grouped_cov_ratios.bed \
  ${out_dir}/${prefix}.merge_comps2.${context}.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.less.grouped.bed_to_join \
  ${out_dir}/${prefix}.use_for_parsing.${context}_intersects.greater.grouped.bed_to_join \
  ${out_dir}/${prefix}.${context}.final_parse \
  ${out_dir}/${prefix}.${context}.final_parse.1.left \
  ${out_dir}/${prefix}.${context}.final_parse.1.right \
  ${out_dir}/${prefix}.${context}.final_parse.1.left.intersect \
  ${out_dir}/${prefix}.${context}.final_parse.1.right.intersect \
  ${out_dir}/${prefix}.merge_comps3.${context}.bed \
  ${out_dir}/${prefix}.merge_comps4.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.${context}_intersect.greater.coord.reduced.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.expanded.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.expanded.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.filtered.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.filtered.${context}.bed \
  ${out_dir}/${prefix}.complement.L_IF_2kb.to_merge.${context}.bed \
  ${out_dir}/${prefix}.complement.R_IF_2kb.to_merge.${context}.bed \
  ${prefix}.allc.${context}_less0.2.bed

  for IF_2kb in L_IF_2kb R_IF_2kb; do
    rm \
    ${out_dir}/${prefix}.complement.${IF_2kb}.filtered.intersect.${context}.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.${context}.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.less.grouped.bed.to_join \
    ${out_dir}/${prefix}.complement.${IF_2kb}.mixed.intersect_total_${context}.greater.grouped.bed.to_join \
    ${out_dir}/${prefix}.complement.${IF_2kb}.to_parse.${context}.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}_intersect.greater.coord.reduced.bed.intersect_less_0.2.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.expanded.${context}_intersect.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.bed \
    ${out_dir}/${prefix}.complement.${IF_2kb}.intersect_total_${context}.grouped_cov_ratios.bed
  done

done

#############################################################################################################################
#############################################################################################################################

# The loop is now complete: integrate the chg and cgn UMRs

# concatenate and merge the chg and chgcgn UMRs
cat ${out_dir}/${prefix}.UMRs_final_merge.chg.bed ${out_dir}/${prefix}.UMRs_final_merge.chgcgn.bed \
| sort -b -k1,1 -k2,2n -k3,3n --parallel=${PBS_NUM_PPN} \
| bedtools merge -d 5 -i stdin > ${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.bed

# Now do the final intersects to get methylation info
allc_chg=${allc_bed_home}/${prefix}.allc.bed.chg

bedtools map -a ${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.bed -b ${allc_chg} -c 6,7,1 -o sum,sum,count -null "0" \
> ${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.mapped.bed

# Now filter and split by size
# removed (($4/$5) < 0.2) && ($5 >= 20) && ($6 >= 5) &&

awk '{ if ( \
($1 != "Pt") && ($1 != "Mt") && ($5 > 0) \
) {print $0,(100*($4/$5))} \
else if ( \
($1 != "Pt") && ($1 != "Mt") && ($5 == 0) \
) {print $0,"-1"}}' OFS="\t" \
${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.mapped.bed > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.total.bed

awk '{ if ( \
(($3-$2) > 150) && ($1 != "Pt") && ($1 != "Mt") && ($5 > 0) \
) {print $0,(100*($4/$5))} \
else if ( \
  (($3-$2) > 150) && ($1 != "Pt") && ($1 != "Mt") && ($5 == 0) \
  ) {print $0,"-1"}}' OFS="\t" \
${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.mapped.bed > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.above150.bed

awk '{ if ( \
(($3-$2) <= 150) && (($3-$2) > 50) && ($1 != "Pt") && ($1 != "Mt") && ($5 > 0) \
) {print $0,(100*($4/$5))} \
else if ( \
  (($3-$2) <= 150) && (($3-$2) > 50) && ($1 != "Pt") && ($1 != "Mt") && ($5 == 0) \
  ) {print $0,"-1"}}' OFS="\t" \
${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.mapped.bed > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.below150.bed

#############################################################################################################################
#############################################################################################################################

# Add info for if UMR is genic, proximal (within 5 kb), or distal (> 5 kb from genes)
for size in above150 below150 total; do

  bedtools intersect -wa -a ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.bed -b ${contigs_no_genes} \
  | awk '{print $0,"-1","distal"}' OFS="\t" \
  > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_no_genes.bed

  bedtools closest -d -t first \
  -a ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.bed \
  -b ${total_genes} \
  | awk '{ if ($14 == 0) {print $1,$2,$3,$4,$5,$6,$7,$14,"genic" } \
  else if ( ($14 > 0) && ($14 <= 5000) ) {print $1,$2,$3,$4,$5,$6,$7,$14,"proximal" } \
  else if ($14 > 5000) {print $1,$2,$3,$4,$5,$6,$7,$14,"distal" } }' OFS="\t" \
    > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_with_genes.bed

  cat \
  ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_with_genes.bed \
  ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_no_genes.bed \
  | sort -b -k1,1 -k2,2n -k3,3n \
  > ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.genic-prox-dist.${size}.bed

done

# Remove temp files

#rm ${out_dir}/${prefix}.UMRs_final_merge.chg.bed
#rm ${out_dir}/${prefix}.UMRs_final_merge.chgcgn.bed
rm ${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.bed
#rm ${out_dir}/${prefix}.UMRs.catted.merged.chg.chgcgn.mapped.bed

for size in above150 below150 total; do
  rm ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_no_genes.bed
  rm ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.${size}.on_contigs_with_genes.bed
done
