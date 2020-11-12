#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
source /work/LAS/mhufford-lab/arnstrm/programs/sourceme
#cp /work/LAS/mhufford-lab/arnstrm/programs/gm_key_64 ~/.gm_key
conda activate braker
genome="$1"
bam="$2" # hintsfile.gff
proteins="$3" #protein_alignment_gth.gff3
#profile="$3"
nam=$(echo ${bam%.*} | sed 's/merged_//g')
today=$(date +"%Y%m%d")
profile=${nam}_prot-rna_${today}.2
GENEMARK_PATH="/work/LAS/mhufford-lab/arnstrm/programs/genemark-et-4.33/bin"
braker.pl \
   --genome=${genome} \
   --hints=${bam} \
   --prot_aln=${proteins} \
   --species=${profile} \
   --softmasking \
   --cores 36 \
   --gff3

