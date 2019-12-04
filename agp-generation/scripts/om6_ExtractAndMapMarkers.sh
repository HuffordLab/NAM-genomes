#!/bin/bash
echo "you provided $# arguments"
if [ "$#" -ne 3 ] ; then
echo "please provide:"
echo -e "\t\t(1) mapped marker file in csv format"
echo -e "\t\t(2) genome in fasta format"
echo -e "\t\t(3) NAM line name"
echo "";
echo "./05_process_runALLMAPS.sh <csv_file> <fasta_file> <NAM-name>" ;
echo "";
exit 0;
fi
module use /work/GIF/software/modules
module load GIF/miniconda/4.5.0
source activate allmaps
module load GIF2/UCSC
module load GIF/allmaps
csv="$1"
fasta="$2"
nam="$3"
# generate pseudo-molecules
python -m jcvi.assembly.allmaps merge ${csv} -o ${csv%.*}.bed
python -m jcvi.assembly.allmaps path ${csv%.*}.bed $fasta
# estimate gap sizes
#python -m jcvi.assembly.allmaps estimategaps ${csv%.*}.bed
#mv ${csv%.*}.estimategaps.agp ${csv%.*}.agp
#python -m jcvi.assembly.allmaps build ${csv%.*}.bed ${fasta}

