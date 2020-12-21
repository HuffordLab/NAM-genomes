#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=30:00:00
#SBATCH --job-name=SNPphyloO1
#SBATCH --output=SNPphyloO1_%j.out 
#SBATCH --error=SNPphyloO1_%j.err
#SBATCH --mail-user=davidehuf@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

ulimit -s unlimited
module load parallel
module load snphylo
cd /work/LAS/mhufford-lab/davehuf/NAMfoundersPrep/SNPdata/Q_SNPhyloKI3/1_rep1

xvfb-run /opt/rit/app/snphylo/20160204/bin/snphylo.sh -H ki3SNPs_merged.hmp -l 0.5 -m 0.05 -M 0.1 -o TIL11 -A -b

