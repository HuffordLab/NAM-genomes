#!/bin/bash

module load blast-plus

# have all your masked blast databases within the same directory

for database in *.nhr

        do

blastn -task dc-megablast -outfmt 6 -query Sbicolor_313_v3.1_exons_primary_notandems_cshl_clusters4.fa -db ${database%.*} -num_threads 4 -out "${database}_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt"

        done
