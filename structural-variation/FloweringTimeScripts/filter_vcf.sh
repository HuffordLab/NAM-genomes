#!/bin/bash

module load vcftools
input=$1
bed=$2

vcftools --vcf ${input} --out ${bed%_*}_SV --bed ${bed} 

#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue

bash scripts/filter_vcf.sh NAM-structural-variations-v2.0.vcf V5_FT_Coordinates.bed
#bash scripts/filter_vcf.sh NAM-structural-variations-v2.0.vcf promoter_V5_FT_Coordinates.bed

#Had to do:
module load bedtools2
tail -n +2 nochr_promoter_V5_FT_Coordinates.bed > tmp.bed
bedtools intersect -header -a tmp.bed -b nochr_NAM-structural-variations-v2.0.vcf  -wa -wb  > intersect_promoter_SVs.bed
tail -n +2 nochr_V5_FT_Coordinates.bed > tmp.bed
bedtools intersect -header -a tmp.bed -b nochr_NAM-structural-variations-v2.0.vcf  -wa -wb  > intersect_V5_SVs.bed
rm tmp.bed