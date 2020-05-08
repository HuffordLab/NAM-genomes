gff2bed < Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff > b73.bed
module use /work/GIF/software/modules
ml GIF/datamash
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose |sed 's/\tN\/A//g' > new-exon-length.txt
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$NF"\t"$0}'   |cut -f 1-4,6- > new-exon-length-with-avg.txt
