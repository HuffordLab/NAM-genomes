#!/bin/bash
genome=$1
module load gmap-gsnap
nam=$(basename $genome |cut -f 1 -d ".")
mkdir -p /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado
gmap_build -d $nam -D /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado ${genome}
gmap \
   -D /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado \
   -d ${nam} \
   -B 4 \
   -t 36 \
   -f gff3_match_cdna \
      /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Trinity/trinity_out_dir/Trinity-GG.fasta \
    > /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_TrinityGG-mapped.gff3  \
    2> /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/stderr-out
