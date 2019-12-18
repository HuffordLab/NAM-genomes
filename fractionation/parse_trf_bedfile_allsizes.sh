#!/bin/bash

for sample in *.fasta.masked.2.7.7.80.10.50.2000.dat
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked.2.7.7.80.10.50.2000.dat//')
                echo $describer

#Convert trf .dat file to .bed file using trf's TRFdat_to_bed.py script:

python TRFdat_to_bed.py --dat ${sample} --bed ${describer}.bed

#Parse above bed file:

cut -f 1,2,3 ${describer}.bed > ${describer}_trf_allsizes.bed

done
