#!/bin/bash
nam=$1
cd /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/$nam/Wallace_chr;
runOM.sh $nam /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_chr/${nam} B73x*_cleaned.vcf /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_chr/${nam}.CHROMOSOMES-broken.fa
cd /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/GoldenGate_chr
runGG.sh $nam /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_chr/${nam}  /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_chr/${nam}.CHROMOSOMES-broken.fa
mkdir -p /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/combined_chr
cp /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/GoldenGate_chr/*_GG-mapped.csv /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/combined_chr/
cp /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/Wallace_chr/*_mapped_onemap-cM.csv /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/combined_chr/
cd /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/combined_chr/
singularity shell --bind $PWD /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/jcvi.simg /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/98_runALLMAPS.sh  *_GG-mapped.csv *_mapped_onemap-cM.csv /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_chr/${nam}.CHROMOSOMES-broken.fa

