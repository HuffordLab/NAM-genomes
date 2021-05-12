#!/bin/bash
#SBATCH --job-name=macs2_call_acrs.nam_lines.v5
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80gb

# macs2_call_acrs.nam_lines.v5.sh

ml MACS2/2.2.7.1-foss-2019b-Python-3.7.4
ml BEDTools/2.29.2-GCC-8.3.0

out_dir=output_directory
tbam_home_dir=directory_containing_bam_files
mapping_list=genome_mapping_combinations.txt
genome_parent_directory=genome_parent_directory
UMR_home_dir=directory_containing_UMR_bed_files_from_step_8
QVALUE=0.005
GENOME_SIZE=1.8e+9

cd ${out_dir}

while read LINE; do

  INPUT_GENOME=$(echo ${LINE} | cut -d ' ' -f1)
  RECIPIENT_GENOME=$(echo ${LINE} | cut -d ' ' -f2)
  COGNATE=$(echo ${LINE} | cut -d ' ' -f3)
  contigs_no_genes=${genome_parent_directory}/${RECIPIENT_GENOME}/processed_genome/${RECIPIENT_GENOME}.chromosomes_contigs_no_genes.bed
  blacklist=${genome_parent_directory}/${RECIPIENT_GENOME}/processed_genome/${RECIPIENT_GENOME}.reduced_header.N_gaps.bed
  total_genes=${genome_parent_directory}/${RECIPIENT_GENOME}/processed_genome/${RECIPIENT_GENOME}.genes.bed

  if [[ $COGNATE == "cognate" ]] ; then
    UMR=${UMR_home_dir}/meth_${INPUT_GENOME}.ref_${INPUT_GENOME}.UMR.above150.with_proximity.bed
  elif [[ $COGNATE == "non-cognate" ]] ; then
    UMR=${UMR_home_dir}/meth_${INPUT_GENOME}.ref_B73.UMR.above150.with_proximity.bed
  else
    echo "error"
  fi

  macs2 callpeak --format BAMPE \
  --gsize ${GENOME_SIZE} \
  -t ${tbam_home_dir}/${INPUT_GENOME}.rep*.${RECIPIENT_GENOME}.dups_removed.mapq30.bam \
  --keep-dup all \
  --qvalue ${QVALUE} \
  --outdir ${out_dir} \
  --name ${INPUT_GENOME}.${RECIPIENT_GENOME}

  # remove blacklist overlaps
  bedtools intersect -wa -v \
  -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.narrowPeak \
  -b ${blacklist} \
  > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak

  # Add info for if ACR is genic, proximal (within 5 kb), or distal (> 5 kb from genes)
  bedtools intersect -wa -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak -b ${contigs_no_genes} \
  | awk '{print $0,"-1","distal"}' OFS="\t" \
  > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.on_contigs_no_genes.bed

  bedtools closest -d -t first \
  -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak \
  -b ${total_genes} \
  | awk '{ if ($17 == 0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17,"genic" } \
  else if ( ($17 > 0) && ($17 <= 5000) ) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17,"proximal" } \
  else if ($17 > 5000) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$17,"distal" } }' OFS="\t" \
  > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.on_contigs_with_genes.bed

  cat \
  ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.on_contigs_with_genes.bed \
  ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.on_contigs_no_genes.bed \
  | sort -b -k1,1 -k2,2n -k3,3n \
  > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.gene_proximity.narrowpeak

  # Calculate intersects with UMRs
  for PROX in genic proximal distal; do
    total_count=`awk -v PROX="$PROX" '{ if ($12==PROX) {print $0} }' OFS="\t" ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.gene_proximity.narrowpeak | sort -u | wc -l`
    overlap_count=`awk -v PROX="$PROX" '{ if ($12==PROX) {print $0} }' OFS="\t" ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.gene_proximity.narrowpeak | bedtools intersect -wa -a stdin -b ${UMR} | sort -u | wc -l`
    echo -e "${overlap_count}\t${total_count}\t${PROX}\t${INPUT_GENOME}\t${COGNATE}" >> ${out_dir}/prox_overlap_rates.temp.txt
  done

done < ${mapping_list}

# calculate the percent of ACR-UMR overlap
awk '{ if ($2 > 0) {print (100*($1/$2)),$0} else {print "-1",$0} }' OFS="\t" ${out_dir}/prox_overlap_rates.temp.txt > ${out_dir}/prox_overlap_rates.final.txt
