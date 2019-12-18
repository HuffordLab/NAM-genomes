#!/bin/bash

module load blast-plus

#looping blastdb:

for sample in *.fasta.masked
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked//')
                echo $describer

makeblastdb -in ${sample} -dbtype nucl -out ${describer}_masked

done
