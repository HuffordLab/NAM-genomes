#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=alignment
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=40:00:00
#SBATCH --mem=40gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jl03308@uga.edu

#files=$(ls /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.fq.gz |cut -f1 -d ".")
genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/B73.PLATINUM.pseudomolecules-v1.fasta
for file in /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.subsample.fq.gz
do
  prefix=$(basename $file | cut -f1 -d ".")
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/SE_alignment.sh $file $genome
done
