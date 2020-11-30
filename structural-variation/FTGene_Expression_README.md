#README for expression data analysis

#Getting gene ids from NAM orthologous to the gene ID's in B73v5

## find NAM coordinates for B73 genes

0. get all gene fastas into the same file

Many of the transcript fastas are only reading as 1 line despite being 2 lines in vim
Must manually go in and add a return; dos2unix won't work

```
cat *.fasta >> FT_candgenes.fasta
```

1. blast the B73 sequences against the NAM genomes
To make all the commands to run the blast jobs concurrently
```
for i in ref_genomes/*.fasta ; do
if [[ ${i} != *B73* ]] ; 
then
	echo bash runblastn.sh ${i} transcript_fastas/FT_candgenes.fasta >> commands.txt ; 
	fi;
done
```
Get 25 commands (1 for each genome of blasting)
Make all the slurm scripts with 
```
python scripts/makeSLURMs.py 1 commands.txt
```
Got multiple hits back, modifying blastn to have a `-max_hsps 1` 
```
cut -f 1 out_B97.pseudomolecules-v1_FT_candgenes | uniq > transcript_name_ref.txt 
#Trying to get just the top hit
awk 'FNR == NR { name[$1] = 0; }
     FNR != NR { for (i in name) if ($0 ~ i && name[i]++ == 0) { print $0; break; } }' \
    reference.txt file.txt
for i in out*candgenes ; do bash gettophits.sh transcript_name_ref.txt ${i} ;done
```
2. get the NAM coordinates for those sequences
Need them in bed file format for bedtools intersect
```
#chr	start	stop	ID
awk -v OFS='\t' '{print $2,$7,$8,$1}' tophitsfile >> bedfile
```
Strandedness is giving me problems because it's not recorded.
```
awk -v OFS='\t' '{if($2 > $3) print $1,$3,$2,$4 ; else print $0}' tophitbedfile 
for i in *tophit_out_*_coords.bed ; do
	awk -v OFS='\t' '{if($2 > $3) print $1,$3,$2,$4 ; else print $0}' ${i} > ${i%coords.bed}nostrd.bed ;
done
```

## find associated NAM gene ID for those coordinates
GFF Files are in `/ptmp/LAS/arnstrm/tpm-final/GFF-PASA` on Condo (released versions)
```
#create soft links
for i in /ptmp/LAS/arnstrm/tpm-final/GFF-PASA/*.gff ; do ln -s ${i} ${i#/ptmp/LAS/arnstrm/tpm-final/GFF-PASA/} ; done

tophit=$1
genome=${tophit#tophit_out_}
genome=${genome%_coords.bed}

GFF=$(ls *.gff | grep ${genome})
module load bedtools2
bedtools intersect -wa -wb -header -a ${GFF} -b ${tophit}

bash filterGFF.sh tophit_out_B97_nostrd.bed

grep "ID=gene" intersect_B97.out |cut -f 10-13 | sort -k1,1 | uniq -c

input=$1
grep "ID=gene" ${input} | cut -f 9 | cut -f 1 -d \; | cut -f 2 -d : > namgeneid.txt
grep "ID=gene" ${input} | cut -f 13 > B73geneid.txt
paste B73geneid.txt namgeneid.txt > ${input%.out}_geneid.txt

bash getgeneid.sh intersect_B97.out
```

## create a matrix for gene IDs
`B73_geneID	NAM1	NAM2	NAM3...`
```
sort -k1,1 intersect_HP301_geneid.txt > test_HP301.txt
sort -k1,1 intersect_P39_geneid.txt > test_P39.txt

#trying Awk
echo B73 HP301 > temp.txt
awk ' 
{ key = $1 } #values it's looking for
!seen[key]++ { keys[++total] = key }
{ values[key] = ( key in values ? values[key] FS $2 : $2 ) }
END {
    for (cnt=1; cnt<=total; cnt++) 
    print keys[cnt], values[keys[cnt]] }' test_HP301.txt > temp.txt
cut -f 1 -d " " temp.txt > temp1.txt
cut -d " " -f 2- temp.txt | awk -v OFS=":" '{i=$(NF); print $1,$2,$3}' > temp2.txt
paste temp1.txt temp2.txt > collapsed.txt

#joining awked collapsed files
join -j 1 -a 1 -a 2 -o auto --header -e NA collapsed_CML228.txt collapsed_B97.txt

join -j 1 -a 1 -a 2 -o auto --header -e NA collapsed_B97.txt collapsed_CML103.txt > temp.txt

join -j 1 -a 1 -a 2 -o auto --header -e NA temp.txt collapsed_CML228.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML247.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML277.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML322.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML333.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML52.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_CML69.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_HP301.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Il14H.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Ki11.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Ki3.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Ky21.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_M162W.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_M37W.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Mo18W.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Ms71.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_NC350.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_NC358.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Oh43.txt |
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Oh7B.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_P39.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Tx303.txt | 
join -j 1 -a 1 -a 2 -o auto --header -e NA - collapsed_Tzi8.txt > allgeneids.txt
```
#Getting the TPM values with the different gene IDs

##Make tissue specific TPM files

Need bam files from LSS for the floral tissues
```
#rsync to /lss/research/mhufford-lab/arnstrm/RNAseq-BAM-NAM-final_2020-03-20
# *V18_ear_*, *V18_tassel*, *R1_anther*
rsync -av /lss/research/mhufford-lab/arnstrm/RNAseq-BAM-NAM-final_2020-03-20/*V18_ear_* /work/LAS/mhufford-lab/snodgras/NAM_FT_genes/RNAseq_bam_files/

#Rename IL14H, KI11, TX303, TZI-8, CML227, Mo18w,Oh7b to Il14H, Ki11, Tx303, Tzi8, CML277, Mo18W, Oh7B

```
```
module load subread
featureCounts -a GFF-FILE -o OUTPUT-counts.txt -T 16 -O -g Parent bamfile1 bamfile2 bamfile3....
```
This was written out in `mkcountfiles.sh`
Then use `calculate-tpm-rpkm-from-feature-counts.R` to create the TPM files
```
ml r-devtools r-dplyr r-tidyr
./calculate-tpm-rpkm-from-feature-counts.R RNAseq_bam_files/B97-counts.txt
```

##Filter the TPM files for NAM specific gene IDs
Run `getTPM.sh` to create the TPM files with only the candidate genes
Then run `formatTPM.sh` to turn the TPMs into long format (detailed below)
```
for n in {1..12} ; 
do let z=19-$n ; 
awk -v OFS='\t' -v tis=$(head -n 1 B97_candTPM.txt | cut -f $z) -v tpm=$(echo ${n}) 
'{print $1,$2,$3,$4,tis,$(NF - tpm)}' B97_candTPM.txt > temp.txt; tail -n +2 temp.txt >> finaltest.txt ; done

input=$1
a=$(awk '{print NF}' $input | head -n 1) #number of columns
head -n 1 $input | cut -f 1-4 > header ; 
awk -v OFS='\t' '{print $0,"Tissue_ID","TPM"}' header > formatted_${input} ; 
for n in $(eval echo "{0..$a}") ; do 
	let z=a-n ; 
	if (( ${z} > 6 )) ; then 
	awk -v OFS='\t' -v tis=$(head -n 1 ${input} | cut -f ${z}) -v tpm=$(echo ${n}) '{print $1,$2,$3,$4,tis,$(NF-tpm)}' $input > temp2.txt ;  
	tail -n +2 temp2.txt >> formatted_${input} ;
	fi ;
done
rm header temp2.txt 

```
Then concatenate those long format files into `allTPM.txt` #need to rm headers so not repeated
```
head -n 1 formatted_B73_candTPM.txt > allTPM.txt
for i in formatted_*; do tail -n +2 $i >> allTPM.txt ; done
```

Add columns for Genome, Tissue, and Replicate by dividing the 5th field, `Tissue_ID`
```
echo Genome > genome.txt
echo Tissue > tissue.txt
echo Replicate > replicate.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 2 -d _ >> genome.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 3,4 -d _ >> tissue.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 5 -d _ >> replicate.txt
sed -i 's/.txt//g' replicate.txt

paste allTPM.txt genome.txt tissue.txt replicate.txt >> allTPM_splitdescriptions.txt
```

### Change Il14H to IL14H and Ms71 to MS71
```
#sed 's/Il14H/IL14H/g' allTPM.txt > allTPM.tmp1
#sed 's/Ms71/MS71/g' allTPM.tmp1 > allTPM.tmp2
#mv allTPM.txt  allTPM.wrongcapitals.txt
#mv allTPM.tmp2 allTPM.txt
#rm allTPM.tmp1 
```


##Create a file with all the TPM values, keep B73v3 GeneIDs

##For the case studies, make sure that all the loci are from the same chromosome

##Make the boxplots with the haplotype vs TMM and then add on the points with jitter

##Ask Ruben what stats he does with his phospholipidase SV study (stat distributions with differing sample sizes)

## Creating the Read pile up version of the TPM
```
module load bedtools2
cat GFF | awk -v OFS="\t" '{if($3 == "gene")print $1,$4,$5,$9}' | bedtools sort | bedtools coverage -a - -b BAM -wao > OUT

for i in B73 B97 ... Tzi8 ; do 
	for j in V11_base V11_middle V11_tip V18_tassel V18_ear R1_anther ; do 
		echo Chr	Start	Stop	ID	counts_${i}_${j}*.bam > ${i}_${j}*.bam_counts.txt
		cat Zm-${i}*.gff | awk -v OFS="\t" '{if($3 == "gene")print $1,$4,$5,$9}' | bedtools sort | bedtools coverage -a - -b ${i}_${j}*.bam -wao > temp.txt
		cut -f 1-5 temp.txt >> ${i}_${j}*.bam_counts.txt
	done
done

cat Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff | awk -v OFS="\t" '{if($3 == "gene")print $1,$4,$5,$9}' | bedtools sort | bedtools coverage -a - -b tissue_bamfiles/CML228_R1_anther_MN06081_R1_round-2Aligned.sortedByCoord.out.bam -wao > temp.txt
cut -f 1-5 temp.txt >> cnts_CML228_R1_anther_MN06081_R1_.txt
rm temp.txt
```
Have to specify each of the bam files individually because there's no good way to do it with the id tag for replicates
This is detailed out in the `mkcountfiles.sh` but is essentially what is outlined above. 

Joining the files using Arun's commons script: `https://github.com/ISUgenomics/common_scripts/blob/master/join_files.sh`
```
sed -e 's/ /\t/g' cnts_B73_R1_anther_MN01081_R1_.txt | awk -v OFS='\t' '{print $4"\;length="$3-$2,$5}' - | tail -n +2 > test1.txt

for i in B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 Oh43 Oh7B P39 Tx303 Tzi8  ; do
	mkdir ${i}_cnts
	for j in raw_bed_cnts/cnts_${i}* ; do
		nam=${j#raw_bed_cnts/}
		sed -e 's/ /\t/g' $j | awk -v OFS='\t' '{print $4"\;length="$3-$2,$5}' - | tail -n +2 > ${i}_cnts/${nam%_R1_.txt}.txt 
	done
	cd  ${i}_cnts/
	bash ../join_files.sh *.txt >> ${i}_joined_cnts.txt
	cd ..
	echo done joining $i
done
```
CML228 joining is having an issue with `CML228_cnts/cnts_CML228_V18_tassel_MN06062.txt`
It has values, but does not have those values in the joined file.
The problem child is  `cnts_CML228_R1_anther_MN06081.txt` which is empty
```
cat Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff | awk -v OFS='\t' '{if($3 == "gene")print $1,$4,$5,$9}' | bedtools sort | bedtools coverage -a - -b tissue_bamfiles/CML228_R1_anther_MN06081_R1_round-2Aligned.sortedByCoord.out.bam  -wao > temp.txt ; cut -f 1-5 temp.txt >> cnts_CML228_R1_anther_MN06081_R1_.txt
```
 

Taking the length information and adding it as a separate column
```
for i in B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 Oh43 Oh7B P39 Tx303 Tzi8  ; do
	cd ${i}_cnts/
	for j in *joined*.txt ; do 
		echo Length > length.txt
		cut -f 1 ${j} | cut -f 4 -d \; | cut -f 2 -d = | tail -n +2 >> length.txt
		paste length.txt ${j} > ${j%.txt}_wLengths.txt
		sed -i 's/\t$//g' ${j%.txt}_wLengths.txt
	done
	cd ..
done
```
Manually add "samples" to the empty 2nd field of the 1st row.

Changed the calculate...R script so it starts at column 2 instead of 7
```
ml r-devtools r-dplyr r-tidyr
./calculate-tpm-rpkm-from-feature-counts.R B73_cnts/B73_joined_cnts.txt
```
The above will not work on a free compute node.

Filter the TPM files for NAM specific gene IDs

Run `getTPM.sh` to create the TPM files with only the candidate genes
```
for i in B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 Oh43 Oh7B P39 Tx303 Tzi8  ; do
	bash /work/LAS/mhufford-lab/snodgras/NAM_FT_genes/scripts/getTPM.sh ${i}_joined_cnts_wLengths_tpm.txt
done
```
Then run `formatTPM.sh` to turn the TPMs into long format
```
input=$1
a=$(awk '{print NF}' $input | head -n 1) #number of columns
head -n 1 $input | cut -f 1-2 > header ; 
awk -v OFS='\t' '{print $0,"Tissue_ID","TPM"}' header > formatted_${input} ; 
for n in $(eval echo "{0..$a}") ; do 
	let z=a-n ; 
	if (( ${z} > 2 )) ; then 
	awk -v OFS='\t' -v tis=$(head -n 1 ${input} | cut -f ${z}) -v tpm=$(echo ${n}) '{print $1,$2,tis,$(NF-tpm)}' $input > temp2.txt ;  
	tail -n +2 temp2.txt >> formatted_${input} ;
	fi ;
done
rm header temp2.txt

for i in *candTPM.txt ; do bash formatTPM.sh $i ;done
```
Then concatenate those long format files into `allTPM.txt` #need to rm headers so not repeated
```
head -n 1 formatted_B73_candTPM.txt > allTPM.txt
for i in formatted_*; do tail -n +2 $i >> allTPM.txt ; done
```

Add columns for Genome, Tissue, and Replicate by dividing the 3rd field, `Tissue_ID`
```
echo Genome > genome.txt
echo Tissue > tissue.txt
echo Replicate > replicate.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 2 -d _ >> genome.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 3,4 -d _ >> tissue.txt
tail -n +2 allTPM.txt | cut -f 3 | cut -f 5 -d _ >> replicate.txt
sed -i 's/.txt//g' replicate.txt

paste allTPM.txt genome.txt tissue.txt replicate.txt >> allTPM_splitdescriptions.txt
``