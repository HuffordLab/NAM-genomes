#!/bin/bash

# extract_NAM_DMRs.v1.sh

# This script finds the inter-NAM DMRs using the conservative split quartile approach
# This script only identifies the DMRs but does not process the gene expression data.

# Explanation of this script

# 1. The coordinates of UMRs in B73 are split into quarters
# 2. The mCHG data from allc files in each NAM line is mapped to the B73 UMR quartiles
# 3. Quartiles meeting the mapping criteria are kept
# 4. B73 UMRs with all four quartiles below 20% mCHG or above 50% mCHG are kept, and designated conserved hypomethylated, or DMR, respectively

out_dir= # directory for output files
merged_list=NAM_mapping_list.txt
allc_chg_home= # directory containing allc methylome files for the CHG context only. Example file: meth_Ki3.ref_B73.allc.bed.chg. File generated in step 6
UMR_home=directory_with_UMR_bed_files_created_in_step_8

# Descriptions of input files:
# ${UMR_home}/${PRIMARY}.UMR.above150.with_proximity.bed
#   example file: meth_B73.ref_B73.UMR.above150.with_proximity.bed
#   col1: chrom, col2: UMR start, col3: UMR end, col4: UMR_name, col5: percent mCHG per UMR, col6: ".", col7: "genic", "proximal", or "distal"
# ${allc_chg_home}/${SECONDARY}.allc.bed.chg
#   example file: meth_Ki3.ref_B73.allc.bed.chg
#     col1: chrom, col2: nucleotide coordinate start (base 0), col3: nucleotide coordinate end (base 1), col4: strand, col5: "CHG", col6: mCHG count, col7: total read count

cd ${out_dir}

while read PRIMARY; do

  # Dissect the B73 leaf UMRs into quarters

  awk '{ print $1,$2,$2+(int(($3-$2+1)/4)),$1,$2,$3,$4,$5,$6,$7}' OFS="\t" ${UMR_home}/${PRIMARY}.UMR.above150.with_proximity.bed > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1.bed
  awk '{ print $1,$2+(int(($3-$2+1)/4)),$2+(2*(int(($3-$2+1)/4))),$1,$2,$3,$4,$5,$6,$7}' OFS="\t" ${UMR_home}/${PRIMARY}.UMR.above150.with_proximity.bed > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q2.bed
  awk '{ print $1,$2+(2*(int(($3-$2+1)/4))),$2+(3*(int(($3-$2+1)/4))),$1,$2,$3,$4,$5,$6,$7}' OFS="\t" ${UMR_home}/${PRIMARY}.UMR.above150.with_proximity.bed > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q3.bed
  awk '{ print $1,$2+(3*(int(($3-$2+1)/4))),$3,$1,$2,$3,$4,$5,$6,$7}' OFS="\t" ${UMR_home}/${PRIMARY}.UMR.above150.with_proximity.bed > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q4.bed

  # Concatenate the UMR quarters into single file
  cat \
  ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1.bed \
  ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q2.bed \
  ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q3.bed \
  ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q4.bed \
  | sort -b -k1,1 -k2,2n -k3,3n --parallel=${SLURM_TASKS_PER_NODE} \
  > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.bed


  # Read each line from the NAM_mapping_list.txt file. This maps the ${SECONDARY}.allc.bed.chg file on top of the B73 UMRs
  while read SECONDARY; do

    OUT=${out_dir}/call_DMR.${PRIMARY}_${SECONDARY}.sh

    echo "#!/bin/bash" > ${OUT}
    echo "#SBATCH --job-name=DMR_${SECONDARY}" >> ${OUT}
    echo "#SBATCH --partition=batch" >> ${OUT}
    echo "#SBATCH --time=24:00:00" >> ${OUT}
    echo "#SBATCH --nodes=1" >> ${OUT}
    echo "#SBATCH --ntasks-per-node=1" >> ${OUT}
    echo "#SBATCH --mem=32gb" >> ${OUT}
    echo "" >> ${OUT}
    echo "# This script was produced by ${out_dir}/${0}" >> ${OUT}
    echo "# This script was produced on `date`" >> ${OUT}
    echo "# primary UMR NAM line is ${PRIMARY}" >> ${OUT}
    echo "# secondary UMR NAM line is ${SECONDARY}" >> ${OUT}
    echo "" >> ${OUT}
    echo "ml BEDTools/2.29.2-GCC-8.3.0" >> ${OUT}
    echo "" >> ${OUT}
    echo "cd ${out_dir}" >> ${OUT}
    echo "" >> ${OUT}
    echo "SECONDARY_ALLC_CHG=${allc_chg_home}/${SECONDARY}.allc.bed.chg" >> ${OUT}

    echo "bedtools map -c 6,7,1 -o sum,sum,count -a ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.bed -b \${SECONDARY_ALLC_CHG} \\" >> ${OUT}
    echo "| awk '{ if (\$12 > 0) {print \$0,(100*(\$11/\$12)) }}' OFS=\"\t\" \\" >> ${OUT}
    echo "> ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.bed" >> ${OUT}

    echo "bedtools groupby -i ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.bed -g 4,5,6,7,8,9,10 -c 11,12,13,11,12,13,4,14,14  -o sum,sum,sum,min,min,min,count,min,max \\" >> ${OUT}
    echo "> ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.quartiles_grouped.bed" >> ${OUT}

    echo "for MIN_COVERAGE in 10; do" >> ${OUT}
    echo "for MIN_COUNT in 4; do" >> ${OUT}

        echo "awk -v MIN_COVERAGE=\"\$MIN_COVERAGE\" -v MIN_COUNT=\"\$MIN_COUNT\" '{ if ( (\$15 >= 50) && (\$12 >= MIN_COVERAGE) && (\$13 >= MIN_COUNT) && (\$14 == 4) ) {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,(100*\$8/\$9)} }' OFS=\"\t\" ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.quartiles_grouped.bed \\" >> ${OUT}
        echo " > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.quartiles_grouped.filtered.mcov_\${MIN_COVERAGE}.mcnt_\${MIN_COUNT}.hyper.bed" >> ${OUT}

        echo "awk -v MIN_COVERAGE=\"\$MIN_COVERAGE\" -v MIN_COUNT=\"\$MIN_COUNT\" '{ if ( (\$16 < 20) && (\$12 >= MIN_COVERAGE) && (\$13 >= MIN_COUNT) && (\$14 == 4) ) {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,(100*\$8/\$9)} }' OFS=\"\t\" ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.quartiles_grouped.bed \\" >> ${OUT}
        echo " > ${out_dir}/${PRIMARY}.UMR.above150.with_proximity.q1234.${SECONDARY}.CHG_mapped.quartiles_grouped.filtered.mcov_\${MIN_COVERAGE}.mcnt_\${MIN_COUNT}.hypo.bed" >> ${OUT}

    echo "done" >> ${OUT}
    echo "done" >> ${OUT}

    qsub ${OUT}

  done < <(grep '.ref_B73' ${merged_list} | cut -f1 | sort -u | grep -v "${PRIMARY}" )
done < <(grep 'meth_B73.ref_B73' ${merged_list} | cut -f1 | sort -u )
