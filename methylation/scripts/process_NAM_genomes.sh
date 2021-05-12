#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N process_NAM_genomes
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=12:00:00

# the genome variable is submitted via an external script: submit_process_NAM_genomes.sh

module load bioawk/1.0-foss-2016b
module load Biopython/1.68-foss-2016b-Python-3.5.2
module load BEDTools/2.28.0-foss-2018a

home_dir=genome_parent_directory
genome_list=NAM_founder_genome_list.txt
cd ${home_dir}

mkdir -p ${home_dir}/${genome}/processed_genome
cd ${home_dir}/${genome}/processed_genome

#this alters the fasta header names by replacing 'chr1' to just '1'
sed 's/^>chr/>/' ${home_dir}/${genome}/raw/${genome}.fa > ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.fa

#this extracts the chromosome and contig lengths from the whole-genome fasta file
bioawk -c fastx '{print $name,"1",length($seq)}' OFS="\t" ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.fa \
| bedtools sort -i stdin | awk '{print $1,$3}' OFS="\t" > ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.bedtools_genome_file.tsv

# This makes 100 bp windows, sliding 20 bp, genome-wide
bedtools makewindows -g ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.bedtools_genome_file.tsv \
-w 100 -s 20 > ${home_dir}/${genome}/processed_genome/${genome}.w100_s20.bed

#this finds the N gaps in the chromosomes and contigs
python findgaps.py \
${home_dir}/${genome}/processed_genome/${genome}.reduced_header.fa \
> ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.N_gaps.temp
bedtools sort -i ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.N_gaps.temp \
> ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.N_gaps.bed
rm ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.N_gaps.temp

#this finds the N-gap complement, which are the contiguous sequences (contigs)
bedtools complement -i ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.N_gaps.bed \
-g ${home_dir}/${genome}/processed_genome/${genome}.reduced_header.bedtools_genome_file.tsv \
| awk '{print $0,"chr_"$1"_new_contig_"NR}' OFS="\t" > ${home_dir}/${genome}/processed_genome/${genome}.homemade_contigs.bed
