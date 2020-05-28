gff2bed < Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff > b73.bed
less b73.bed 
less -S b73.bed
cut -f 1-3,6 b73.bed |head
cut -f 1-3,6,7 b73.bed |head
cut -f 1-3,6,8 b73.bed |head
cut -f 1-3,6,8,9 b73.bed |head
cut -f 1-3,6,8,10 b73.bed |head
cut -f 1-3,6,8,10 b73.bed |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |sed 's/;/\t/g' |cut -f1-7 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 2 -d ";" |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";" |less -S
#awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";" |awk '{print $3-$2"\t"
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";" |head
#awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $3-$2"\t"$7}' |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $3-$2"\t"$7}' |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |less -S 
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |head
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon/\texon/g' |head
module use /work/GIF/software/modules
ml GIF/datamash
#datamash -s crosstab 2,1 sum 5 <all-sn-values.txt
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon/\texon/g' |datamash -s crosstab 2,1 sum 3 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon/\texon/g' |datamash -s crosstab 1,2 sum 3 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon/\texon/g' |datamash -s crosstab 1,2 sum 3 > full-table.txt
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon/\t/g' |datamash -s crosstab 1,2 sum 3 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 1,2 sum 3 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose |sed 's/\tN\/A//g' |less -S
awk '$8=="exon"' b73.bed |cut -f 1-3,6,8,10 |cut -f 1-2 -d ";"|sed 's/;/\t/g' |awk '{print $7"\t"$3-$2}' |sed 's/Name=//g' |sed 's/.exon./\t/g' |datamash -s crosstab 2,1 sum 3 |sort -k1,1 -n |datamash transpose |sed 's/\tN\/A//g' > new-exon-length.txt
less -S awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.0f", a[i]/NR; printf "\t"};printf "\n"}'
#awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.0f", a[i]/NR; printf "\t"};printf "\n"}'
lt
less -S new-exon-length.txt 
awk 'NF<=3' new-exon-length.txt |less -S
awk 'NF>3' new-exon-length.txt |less -S
awk 'NF>3' new-exon-length.txt |awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.0f", a[i]/NR; printf "\t"};p
#awk 'NF>3' new-exon-length.txt |awk '{for (i=1;i<=NF-1;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.0f", a[i]/NR; printf "\t"};printf "\n"}'
awk 'NF>3' new-exon-length.txt |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=3;i<=NF-1;i++){a[i]+=$i;}} END {for (i=3;i<=NF-1;i++){printf "%.0f", a[i]/NR; printf "\t"};printf "\n"}' |less -S
lt
nano avg.sh
lt
chmod +x avg.sh 
lt
awk 'NF>3' new-exon-length.txt > temp.txt
./avg.sh temp.txt 
./avg.sh temp.txt |less -S
nano avg.sh
awk 'NF>3' new-exon-length.txt | awk '{for (i=3;i<=NF-1;i++){a[i]+=$i;}} {print $0"\t"a[i]/i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=NF-1;i++){a[i]+=$i;}} {print $0"\t"a[i]/i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=3;i<=NF-1;i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=NF-1;i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"i}' |less -S
#awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=NF-1;i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"i}' |less -S
awk '{print NF,NF-1}'  new-exon-length.txt |less -S
awk '{print NF,(NF-1)}'  new-exon-length.txt |less -S
#awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"$i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=3;i<=(NF-1);i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"$i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){a[i]+=$i;}} {print $0"\t"a[i]"\t"$i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){sum[i] += $i;}} {print $0"\t"sum[i]"\t"$i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{for (i=2;i<=(NF-1);i++){sum[i] += $i;}} {print $0"\tsum: "sum[i]"\t"$i}' |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' input.txt
salloc -p freemedium  -n 16 -t 8:00:00 -N 1
lt
qsa
ll
qsa
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-2; print $0"\t"sum}'  |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-2; print $1"\t"sum"\t"$0}'   |less -S
awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-2; print $1"\t"sum"\t"$0}' new-exon-length.txt |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-2; print $1"\t"sum"\t"$0}'   |cut -f 1-2,4- |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"$0}'   |cut -f 1-2,4- |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$0}'   |cut -f 1-2,4- |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$0}'   |cut -f 1-3,5- |less -S
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$0}'   |cut -f 1-3,5- > new-exon-length-with-avg.txt
awk 'NF>3' new-exon-length.txt | awk '{sum = 0; for (i = 3; i <= NF-1; i++) sum += $i; sum /= NF-3; print $1"\t"sum"\t"NF-1"\t"$NF"\t"$0}'   |cut -f 1-4,6- > new-exon-length-with-avg.txt
less -S new-exon-length-with-avg.txt
salloc -p freemedium  -n 16 -t 8:00:00 -N 1
lt
less temp.txt 
history |tail -n 100 > cmds.sh
