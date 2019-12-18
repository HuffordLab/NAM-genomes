#!/bin/bash

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt//')
                echo $describer


#First, do another filter for tandem repeats using Dagchainer's filter_repetitive_matches.pl script:

perl /DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 100000 < ${sample} > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered;

#Run dagchainer (the output will append the extension .aligncoords to each output file); parameters were optimized for maize, and -A 15 takes into consideration the fact that individual exons, not genes, are being processed:

perl /DAGCHAINER/run_DAG_chainer.pl -i ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered -s -I -D 1000000 -g 40000 -A 15

        done
