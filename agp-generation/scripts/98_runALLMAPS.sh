#!/bin/bash
export PATH=$PATH:/shared/software/GIF/programs/UCSC/323/bin
csv1="$1"
csv2="$2"
fasta="$3"
python -m jcvi.assembly.allmaps merge ${csv1} ${csv2} -o combined.bed
#sed -i 's/_GG_mapped 1/_GG_mapped 2/g' weights.txt
python -m jcvi.assembly.allmaps path combined.bed $fasta

#python -m jcvi.assembly.allmaps estimategaps ${csv%.*}.bed
#mv ${csv%.*}.estimategaps.agp ${csv%.*}.agp
#python -m jcvi.assembly.allmaps build ${csv%.*}.bed ${fasta}

