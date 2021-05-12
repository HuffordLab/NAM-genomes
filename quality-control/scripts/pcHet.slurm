#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --job-name=perHet
#SBATCH --output=perHet_%j.out
#SBATCH --error=perHet_%j.err
#SBATCH --mail-user=davidehuf@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

ulimit -s unlimited
module load parallel
cd /work/LAS/mhufford-lab/davehuf/NAMfoundersPrep/SNPdata/SNPs
python ../18_1_countPerPolyLoci.py maizeHapMapV2_B73RefGenV2_201203028_Allchr.hmp.txt.1title maizeHapMapV2_B73RefGenV2_201203028_Allchr.hmp.txt.1title.perHet

