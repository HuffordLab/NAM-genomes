#!/bin/bash

# This script processes the bowtie2-aligned ATAC-seq reads for all of the NAM lines, mapped to the B73 reference genome and their cognate reference genomes

out_dir=output_directory
mapping_list=genome_mapping_combinations.txt
ATAC_home_dir=ATAC_parent_directory

while read LINE; do

  INPUT_GENOME=$(echo ${LINE} | cut -d ' ' -f1)
  RECIPIENT_GENOME=$(echo ${LINE} | cut -d ' ' -f2)
  OUT=${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}.sort_bam.sh

  echo "#!/bin/bash" > ${OUT}
  echo "#SBATCH --job-name=${INPUT_GENOME}.${RECIPIENT_GENOME}.sort_bam" >> ${OUT}
  echo "#SBATCH --partition=batch" >> ${OUT}
  echo "#SBATCH --time=12:00:00" >> ${OUT}
  echo "#SBATCH --nodes=1" >> ${OUT}
  echo "#SBATCH --ntasks-per-node=12" >> ${OUT}
  echo "#SBATCH --mem=80gb" >> ${OUT}

  echo "ml deepTools/3.3.1-intel-2019b-Python-3.7.4" >> ${OUT}
  echo "ml SAMtools/1.10-iccifort-2019.5.281" >> ${OUT}

  echo "cd ${out_dir}" >> ${OUT}
  echo "for SAM in ${ATAC_home_dir}/${INPUT_GENOME}/aligned/${INPUT_GENOME}.rep*.${RECIPIENT_GENOME}.sam; do" >> ${OUT}
  echo "  sam_basename=\`basename --suffix=.sam \"\$SAM\"\`" >> ${OUT}
  echo "  # convert the SAM to sorted BAM" >> ${OUT}

# convert sam to bam
  echo "  samtools view -b -h -S --threads \${SLURM_TASKS_PER_NODE} \\" >> ${OUT}
  echo "  -o ${out_dir}/\${sam_basename}.unsorted.bam \${SAM}" >> ${OUT}

# remove duplicates
  echo "  sambamba markdup \\" >> ${OUT}
  echo "  --remove-duplicates -t \${SLURM_TASKS_PER_NODE} --tmpdir ${out_dir} \\" >> ${OUT}
  echo "  ${out_dir}/\${sam_basename}.unsorted.bam \\" >> ${OUT}
  echo "  ${out_dir}/\${sam_basename}.dups_removed.bam" >> ${OUT}

# filter by MAPQ >= 30
  echo "  sambamba sort \\" >> ${OUT}
  echo "  --memory-limit=76GB \\" >> ${OUT}
  echo "  --tmpdir ${out_dir} \\" >> ${OUT}
  echo "  -o ${out_dir}/\${sam_basename}.dups_removed.mapq30.bam \\" >> ${OUT}
  echo "  -t \${SLURM_TASKS_PER_NODE} \\" >> ${OUT}
  echo "  -F \"mapping_quality >= 30\" \\" >> ${OUT}
  echo "  ${out_dir}/\${sam_basename}.dups_removed.bam" >> ${OUT}

# Create 1 bp resolution bigwig file
  echo "  bamCoverage --outFileFormat bigwig --binSize 1 --numberOfProcessors \${SLURM_TASKS_PER_NODE} \\" >> ${OUT}
  echo "  -b ${out_dir}/\${sam_basename}.dups_removed.mapq30.bam \\" >> ${OUT}
  echo "  -o ${out_dir}/\${sam_basename}.dups_removed.mapq30.bigwig" >> ${OUT}
  echo "done" >> ${OUT}

qsub ${OUT}
done < ${mapping_list}
