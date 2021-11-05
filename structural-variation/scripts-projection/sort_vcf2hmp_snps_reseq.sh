#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=10,mem=120gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp/
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp/
#PBS -V
#PBS -N sort_vcf2hmp_snps_reseq
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

module load parallel

parallel --jobs 10 < ~/projects/sv_nams/scripts/commands_sort_vcf_reseq.txt
parallel --jobs 10 < ~/projects/sv_nams/scripts/commands_vcf2hmp_reseq.txt
