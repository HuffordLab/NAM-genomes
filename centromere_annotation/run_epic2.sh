#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=trimming
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=40:00:00
#SBATCH --mem=40gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jl03308@uga.edu

#files=$(ls /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.fq.gz |cut -f1 -d ".")
for ChIP in /scratch/jl03308/NAM_pancentromere/analysis/peak_call/*/*.ChIP.q20.sorted.bam
do
  prefix=$(basename $ChIP | cut -f1 -d ".")
  Input=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${prefix}/${prefix}.q20.sorted.bam
  wdir=$(dirname $ChIP)
  cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${prefix}.*pseudomolecules-v*.fasta.fai |cut -f1,2 \
  > ${wdir}/${prefix}.chrom.sizes
  chrsize=${wdir}/${prefix}.chrom.sizes
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/epic2.sh $ChIP $Input $chrsize
done
