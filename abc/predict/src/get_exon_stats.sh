#gff="data/ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff"
gff="data/ref/zea_maysb73_core_3_87_1_chr-added.gff"
window="data/ref/Zm-B73-REFERENCE-NAM-5.0_20Mb.bed"
foldfile="data/ref/Zm-B73-REFERENCE-NAM-5.0_FOLD.txt"

awk '$5 > $4' $gff > TMP

gff="TMP"

echo "gene gaps"
awk '$3 ~ /gene/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | bedtools complement -g data/ref/Zm-B73-REFERENCE-NAM-5.0.gbed -i stdin | awk '{print $3 - $2}' > data/ref/Zm_intergenic.txt

echo "cds gaps"
awk '$3 ~ /CDS/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | bedtools complement -g data/ref/Zm-B73-REFERENCE-NAM-5.0.gbed -i stdin | awk '{print $3 - $2}' > data/ref/Zm_intragenic.txt

echo "cds sizes"
awk '$3 ~ /CDS/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | awk '{print $3 - $2}' > data/ref/Zm_cds_sizes.txt

echo "cds numbers"
awk '$3 ~ /CDS/' $gff | grep -v "scaf" | cut -f9 | sort | uniq -c | awk '{print $1}' > data/ref/Zm_cds_numbers.txt

echo "genes per XMb window"
awk '$3 ~ /gene/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | bedtools intersect -wa -a $window -b stdin | uniq -c | awk '{print $1}' > data/ref/Zm_genes_per_window.txt

echo "bps per XMb window"
awk '$3 ~ /CDS/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | bedtools intersect -wa -wb -a $window -b stdin | awk '{print $0 "\t" $6 - $5}' | bedtools groupby -i stdin -grp 1,2,3 -c 7 -o sum | awk '{print $4}' > data/ref/Zm_bps_per_window.txt

echo "0 and 4 fold counts per window"
awk 'NR > 1 && $3 == 0 {print $1 "\t" $2-1 "\t" $2}' $foldfile | grep -v scaf | bedtools sort -i stdin | bedtools intersect -wb -a stdin -b $window | cut -f4-6 | sort | uniq -c | awk '{print $1}' > data/ref/Zm_0fold_per_window.txt

awk 'NR > 1 && $3 == 4 {print $1 "\t" $2-1 "\t" $2}' $foldfile | grep -v scaf | bedtools sort -i stdin | bedtools intersect -wb -a stdin -b $window | cut -f4-6 | sort | uniq -c | awk '{print $1}' > data/ref/Zm_4fold_per_window.txt

for i in `ls data/ref/Zm_*`; do echo $i; cat $i | Rscript -e 'summary (as.numeric (readLines ("stdin")))'; done


#build data frame with header for files with 113 total lines
header=`wc -l data/ref/Zm_* | awk '$1 == 113 {print $2}' | sed 's/_/\t/g' | cut -f2,3,4 | sed 's/\twindow.txt//g;s/.txt//g;s/\tper//g;s/\t\t*/_/g' | tr "\n" "\t" `

#build data frame
cat <(echo -e "chrom\tstart\t$header") <(paste data/ref/Zm-B73-REFERENCE-NAM-5.0_20Mb.bed `wc -l data/ref/Zm_* | awk '$1 == 113 {print $2}'`) > data/ref/window_stats_20Mb.txt


rm TMP
