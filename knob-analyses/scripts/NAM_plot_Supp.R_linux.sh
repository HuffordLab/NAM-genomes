
###PREP annotation files for supplement figure 8


## plotting was done in local R, so the files in original state are too big
## for this script, all gene annotation files (all ending with 1.gff)  
## and all TE annotation files (all ending with .gff3 ) are in the same directory

module load BEDTools/2.26.0-GCCcore-8.3.0

# TE annotations
ls *.gff3 > fils

while read f
do
	#line to modify col 1
	awk '{print $1}' $f | sed 's/.*_//'  > new_col1

	#line to modify type col
	awk '{print $9}' $f | q | sed 's/;Class.*//'  > new_col2
	awk '{print $9}' $f | sed 's/.*Classification=//' | sed 's/;Seq.*//'  > new_col3

	awk '{ print $4"\t"$5}' $f  > mid
	paste new_col1 mid new_col2 new_col3 > check

	#only the chr line
	sed -n '/^chr/p' check> check_chr.bed

	#get the line name to subset the arrays relavent, then intersect with the annotations
	tar=$(sed -n '7 p' $f | sed 's/_.*//')
	if [[ "$tar" == "Il14H" ]]; then tar="IL14H" ; fi
	if [[ "$tar" == "Ms71" ]]; then tar="MS71" ; fi
	sed -n "/^$tar/p" NAM_array_coords.tsv > sub
	awk '{ print $2"\t"$4"\t"$5}' sub  | sed -n '/^chr/p' > sub.bed
	bedtools intersect -a sub.bed -b check_chr.bed -wb  > reduced_check_chr.bed

	#final file for each line
	awk '{ print $1"\t"$2"\t"$3"\t"$7"_"$8}' reduced_check_chr.bed > mod_$f
done < fils


# gene annotations
ls zea*1.gff > fils2

#convert first column to all upper for easier filtering
awk '$1 = toupper($1)'  NAM_array_coords.tsv > ed_NAM_array_coords.tsv 

#TE files for matching
ls mod*.gff3 > match_TE
awk '$1 = toupper($1)' match_TE > match_TE_col2
paste match_TE  match_TE_col2 > new_match_TE



while read f
do
	awk -F '\t' '$3 == "gene" { print }' $f |  awk '{ print $1"\t"$4"\t"$5"\t"$3}' | sed -n '/^chr/p' > gene.bed
	tar=$(echo $f | sed 's/_core.*//' | sed 's/.*mays//')
	tarup=${tar^^}
	
	sed -n "/^$tarup/p" ed_NAM_array_coords.tsv  | awk '{print $2"\t"$4"\t"$5}'| sed -n '/^chr/p' > sub.bed
	
	bedtools intersect -a sub.bed -b gene.bed -wb  > reduced_gene.bed
	
	awk '{ print $1"\t"$2"\t"$3"\t"$7}' reduced_gene.bed > mod_$f
	
	#now match to the TE file
	match=$(awk "/$tarup/"  new_match_TE | awk '{print $1}')
	
	cat mod_$f $match | bedtools sort -i | uniq > all_$tarup.bed 
	
done < fils2

########


