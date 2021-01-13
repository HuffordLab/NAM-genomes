#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N remove_blacklist.change_filenames.convert_bed.NAM_UMRs.v1
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=1:00:00

# remove_blacklist_NAM_UMRs.v1.sh

# Description:
# (1) Remove the blacklisted regions from the previously called UMRs
# (2) Rename the UMRs to the standardized naming scheme
# (3) Convert the files to the formal BED format.

ml BEDTools/2.28.0-foss-2018b

out_dir=output_directory
blacklist_home_dir=directory_containing_blacklist_bed_files
UMR_home_dir=directory_containing_UMR_files_from_step_7

cd ${out_dir}

mkdir -p ${out_dir}/umrs_with_prox_info

# Now remove blacklist, change to bed format, and convert filenames
while read LINE; do

  prefix=$(echo ${LINE} | cut -d ' ' -f1)
  REF_SHORT_NAME=$(echo ${LINE} | cut -d ' ' -f2)
  REF_LONG_NAME=$(echo ${LINE} | cut -d ' ' -f3)

  bedtools intersect -v -wa \
  -a ${UMR_home_dir}/${prefix}.final_UMRs.chg.chgcgn.genic-prox-dist.above150.bed \
  -b ${blacklist_home_dir}/${REF_SHORT_NAME}_PE_150_30x.l.vsl.off4.seed1.dups_removed.mapq_20.seed_cutoff80.spread_cutoff25.n_gaps_removed.bed \
  | awk -v prefix="$prefix" '{print $1,$2,$3,prefix"_"$1":"$2".."$3,$7,".",$9}' OFS="\t" \
  > ${out_dir}/umrs_with_prox_info/${prefix}.UMR.above150.with_proximity.bed

  bedtools intersect -v -wa \
  -a ${UMR_home_dir}/${prefix}.final_UMRs.chg.chgcgn.genic-prox-dist.total.bed \
  -b ${blacklist_home_dir}/${REF_SHORT_NAME}_PE_150_30x.l.vsl.off4.seed1.dups_removed.mapq_20.seed_cutoff80.spread_cutoff25.n_gaps_removed.bed \
  | awk -v prefix="$prefix" '{print $1,$2,$3,prefix"_"$1":"$2".."$3,$7,".",$9}' OFS="\t" \
  > ${out_dir}/umrs_with_prox_info/${prefix}.UMR.total.with_proximity.bed

  cut -f1-6 ${out_dir}/umrs_with_prox_info/${prefix}.UMR.above150.with_proximity.bed > ${out_dir}/${prefix}.UMR.bed

done < <(cat ${out_dir}/methylomes_table_cross_ref.txt)
