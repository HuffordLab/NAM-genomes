#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=subsampling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem=20gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jl03308@uga.edu

#files=$(ls /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.fq.gz |cut -f1 -d ".")
for file in /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*_trimmed.fq.gz
do
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/SE_subsampling.sh $file 5000000
done
