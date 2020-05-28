#!/bin/bash
ml bedops
gff=$1
gff2bed < $gff > ${gff%.*}.bed
module use /work/GIF/software/modules
ml GIF/datamash
awk '$8=="exon"' ${gff%.*}.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose |sed 's/\tN\/A//g' > new-exon-length.txt
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$NF"\t"$0}'   |cut -f 1-4,6- > new-exon-length-with-avg.txt


for f in B73Ab10.bed Oryza_sativa.bed Sorghum_bicolor.bed ZeamaysV3.bed ZeamaysV4.bed ZeamaysV5.bed zmMAKERPASA.bed; do
awk '$8=="exon"'$f |sed 's/;/\t/g' |cut -f 1-3,6,8,10


for bed in B73Ab10.bed; do
awk '$8=="exon"' $bed |sed 's/;/\t/g' |cut -f 1-3,6,8,10 | sed 's/Parent=//g' | awk '{print $6"\t"$3-$2}' | sed 's/.exon./\t/g' | datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose | sed 's/\tN\/A//g'  > ${f%.*}_exon-stats.txt
awk 'NF>3' ${f%.*}_exon-stats.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$NF"\t"$0}'   |cut -f 1-4,6- >  ${f%.*}_exon-stats-avg.txt
done
for bed in ZeamaysV3.bed; do
awk '$8=="exon"' $bed |sed 's/;/\t/g' |cut -f 1-3,6,8,10 | sed 's/Parent=//g' | awk '{print $6"\t"$3-$2}' | sed 's/_T/\t/2' | datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose | sed 's/\tN\/A//g'  > ${f%.*}_exon-stats.txt
awk 'NF>3' ${f%.*}_exon-stats.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$NF"\t"$0}'   |cut -f 1-4,6- >  ${f%.*}_exon-stats-avg.txt
done



for f in Oryza_sativa.bed Sorghum_bicolor.bed ZeamaysV4.bed ZeamaysV5.bed zmMAKERPASA.bed; do
awk '$8=="exon"' $bed |sed 's/;/\t/g' |cut -f 1-3,6,8,11 |






B73Ab10.bed
chr1    46192   46888   +       exon    Parent=mikado.chr1G1.1
chr1    46192   46888   +       exon    Parent=mikado.chr1G1.2
Oryza_sativa.bed
1       2982    3268    +       exon    Name=Os01t0100100-01-E1
1       3353    3616    +       exon    Name=Os01t0100100-01-E2
Sorghum_bicolor.bed
1       1950    2454    +       exon    Name=EER90453-1
1       2472    2616    +       exon    Name=EER90453-2
ZeamaysV3.bed
1       2       104     +       exon    Parent=GRMZM2G060082_T01
1       199     313     +       exon    Parent=GRMZM2G060082_T01
ZeamaysV4.bed
1       44288   44947   +       exon    Name=Zm00001d027230_T001.exon1
1       45665   45803   +       exon    Name=Zm00001d027230_T001.exon2
ZeamaysV5.bed
chr1    34616   35318   +       exon    Name=Zm00001e000001_T001.exon.1
chr1    36036   36174   +       exon    Name=Zm00001e000001_T001.exon.2
zmMAKERPASA.bed
1       34616   35318   +       exon    Name=Zm00001e000001_T001.exon.1
1       36036   36174   +       exon    Name=Zm00001e000001_T001.exon.2

