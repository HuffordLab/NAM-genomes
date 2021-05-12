#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=trimming
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem=20gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jl03308@uga.edu

cd /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/before_merge/
# lines=$(ls *.gz | cut -f1 -d "_" |sort |uniq)
# for line in $lines
# do
#   cat ${line}_*.gz > ../${line}.fq.gz
# done

for file in /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.fq.gz
do
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/SE_trimming.sh $file
done
