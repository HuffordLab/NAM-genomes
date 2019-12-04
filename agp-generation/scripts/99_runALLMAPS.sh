#!/bin/bash
export PATH=$PATH:/shared/software/GIF/programs/UCSC/323/bin
csv="$1"
fasta="$2"
nam=$(basename $csv |cut -f 1 -d ".")
python -m jcvi.assembly.allmaps merge ${csv} -o ${csv%.*}.bed
python -m jcvi.assembly.allmaps path --mincount=3 --format=png ${csv%.*}.bed $fasta

#python -m jcvi.assembly.allmaps estimategaps ${csv%.*}.bed
#mv ${csv%.*}.estimategaps.agp ${csv%.*}.agp
#python -m jcvi.assembly.allmaps build ${csv%.*}.bed ${fasta}

