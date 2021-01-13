#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N process_gene_annotations
#PBS -l nodes=1:ppn=4
#PBS -l mem=12gb
#PBS -l walltime=12:00:00

module load BEDTools/2.28.0-foss-2018a
genome_list=NAM_founder_genome_list.txt
home_dir=genome_parent_directory

cd ${home_dir}

while read GENOME; do

  out_dir=${home_dir}/${GENOME}/processed_genome

  # convert from gff to bed file, extract the genes, and convert chromosome names
  awk -F"\t|:|;" '{ if ( $3=="gene" ) { print $1,$4,$5,$10,$3,$7} }' OFS="\t" \
  ${home_dir}/${GENOME}/raw/${GENOME}_*.gff3 \
  | sed 's/^chr//' \
  | sort -b -k1,1 -k2,2n -k3,3n  > ${out_dir}/${GENOME}.genes.bed

  genome_file=${home_dir}/${GENOME}/processed_genome/${GENOME}.reduced_header.bedtools_genome_file.tsv

  # convert the genome file to a bed file (temporary)
  awk '{print $1,"1",$2}' OFS="\t" ${genome_file} > ${out_dir}/chromosomes_contigs.temp.bed

  # Find chromosomes and contigs that have no genes on them
  bedtools intersect -v -wa -a ${out_dir}/chromosomes_contigs.temp.bed -b ${out_dir}/${GENOME}.genes.bed \
  | sort -b -k1,1 -k2,2n -k3,3n > ${out_dir}/${GENOME}.chromosomes_contigs_no_genes.bed

  # Find chromosomes and contigs that have genes on them
  bedtools intersect -wa -u -a ${out_dir}/chromosomes_contigs.temp.bed -b ${out_dir}/${GENOME}.genes.bed \
  | sort -b -k1,1 -k2,2n -k3,3n > ${out_dir}/${GENOME}.chromosomes_contigs_with_genes.bed

  # clean up
  rm ${out_dir}/chromosomes_contigs.temp.bed

done < <(cut -f1 ${genome_list} | sort -u)
