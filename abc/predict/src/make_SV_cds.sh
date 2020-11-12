#grep -Ff <(zcat data/variants/NAM-structural-variations-v3.0.txt.gz | cut -f1 | grep -v scaf | sed 's/[-:]/\t/g' | awk 'NR > 1' | awk 'BEGIN{OFS="\t"}{print $0, $1 "-" $2 ":" $3 "-" $4}' | awk 'BEGIN{OFS="\t"}{if($1 != $3){print $1, $2, $1, $2, $5 "\n" $3, $4, $3, $4, $5} else{print $0}}' | awk 'BEGIN{OFS = "\t"}{if($2 < $4){print $1, $2-1, $4, $5}else{print $1,$4-1,$2, $5}}' | bedtools sort -i stdin | bedtools intersect -wa -a stdin -b <(awk '$3 ~ /CDS/' data/ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff | grep -v scaf) | cut -f4 | sort -u) <(zcat data/variants/NAM-structural-variations-v3.0.txt.gz) > data/variants/NAM-structural-variations_CDS.txt

sv="data/variants/NAM_founders_SVs.sniffles-bionano.hmp.txt"
sv_out="data/variants/NAM-structural-variations_CDS.txt"
gff="data/ref/zea_maysb73_core_3_87_1_chr-added.gff"

awk '$5 > $4' $gff > TMP

gff="TMP"



grep -Ff <(bedtools sort -i <(python src/parse_sv.py $sv) | bedtools intersect -wa -a stdin -b <(awk '$3 ~ /CDS/' $gff | grep -v scaf) |  cut -f4 | sort -u) $sv > $sv_out


rm TMP
