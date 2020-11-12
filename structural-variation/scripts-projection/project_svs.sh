#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=20gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/projection
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/projection
#PBS -V
#PBS -N project_svs_${CROSS}_${HAPSIZE}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to data folder
cd ~/projects/sv_nams

# create haplotypes from parents
run_pipeline.pl -Xmx20g -FILLINFindHaplotypesPlugin -hmp data/NAM_parents_SVs-SNPs.${CROSS}.sorted.hmp.txt  -o analysis/projection/donors_${CROSS} -hapSize ${HAPSIZE} -minTaxa 1

# project ril genotypes based on parental haplotypes
run_pipeline.pl -Xmx20g -FILLINImputationPlugin -hmp data/NAM_rils_SVs-SNPs.${CROSS}.best-markers.not-projected.sorted.hmp.txt -d analysis/projection/donors_${CROSS} -o analysis/projection/NAM_rils_SVs-SNPs.${CROSS}.best-markers.projected.hmp.txt -hapSize ${HAPSIZE} -accuracy -hybNN false

# convert projected hapmap to diploid format
run_pipeline.pl -Xmx20g -importGuess analysis/projection/NAM_rils_SVs-SNPs.${CROSS}.best-markers.projected.hmp.txt -export analysis/projection/NAM_rils_SVs-SNPs.${CROSS}.best-markers.projected.hmp.txt -exportType HapmapDiploid
