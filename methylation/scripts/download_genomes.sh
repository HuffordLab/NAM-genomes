#!/bin/bash
# Download the 25 NAM founder genomes and B73 V5

home_dir=genome_parent_directory
genome_list=NAM_founder_genome_list.txt

while read genome; do
  mkdir -p ${home_dir}/${genome}/raw/
  cd ${home_dir}/${genome}/raw/
  wget --recursive --no-parent --no-host-directories --cut-dirs=6 https://download.maizegdb.org/${genome}/ ./
  gunzip ${home_dir}/${genome}/raw/*.gz
done < <(cut -f1 ${genome_list} | sort -u)
