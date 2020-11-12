# Detecting SVs near known flowering time genes

# Data
See `FT_genename_V3id.txt` for gene name from Dong et al 2012 to B73V3 or GenBank accession ID

# Creating SV tracks for candidate regions
## Formatting SV file
```
# replace AA with 0/0(absent) and TT with 1/1(present) and NN with ./.(missing)
sed -e 's+AA+0/0+g' -e 's+TT+1/1+g' -e 's+NN+./.+g' NAM_founders_SVs.sniffles-bionano.hmp.txt > Nam_replaced.txt
# get SV_IDs
cut -f 1 Nam_replaced.txt > SV_id.txt
# get SV_genotypes
cut -f 12- Nam_replaced.txt > Nam_genotypes.txt
# combine two files
paste SV_id.txt Nam_genotypes.txt > combined.txt
# add id number to the beginning of the first field
awk -v OFS='\t' '{print NR-1 $0}' combined.txt > combined_wIDs.txt 
# makeGenome SV files for each genome
bash mkGenomeSV.sh combined_wIDs.txt
# get rid of second column
for i in *combined_wIDs.txt; do cut -f 1 ${i} > tmp_${i} ; done
# get rid of header
for i in tmp*; do tail -n +2 ${i} > tmp2_${i} ; done
# replace . with tab for SV_IDs
for i in tmp2*; do sed 's/\./ /g' ${i} > ${i#tmp2_tmp_}.bed; done
# changes the ending of the file names to .bed
for i in *txt.bed ; do  mv ${i} ${i%.txt.bed}.bed ;done
#remove the temporary files that you don't want and the AB10 files 
rm tmp* B73_Ab10*

#Make individual files for each genome and SV type
for i in *wIDs.bed; do grep ins ${i} > INS_${i} ; done``
mkdir INS
mv INS_* INS/.
for i in *wIDs.bed; do grep del ${i} > DEL_${i} ; done
mkdir DEL
mv DEL_* DEL/.
for i in *wIDs.bed; do grep dup ${i} > DUP_${i} ; done
mkdir DUP
mv DUP_* DUP/.
for i in *wIDs.bed; do grep inv ${i} > INV_${i} ; done
mkdir INV
mv INV_* INV/.

#make them in bed file format (change order of columns)

for j in DEL DUP INV INS ; do cd $j ; 
	for i in *wIDs.bed ; do awk -v OFS='\t' '{print $2,$3,$4,$1}' $i > truebed_${i} ; done; 
	cd .. ; done #get in the correct bed file format
 #clean up directory
for j in DEL DUP INV INS ; do cd $j ; mkdir temp ; mv truebed* temp/. ; rm *wIDs.bed ; mv temp/* . ; rmdir temp ; for i in truebed_* ; do mv $i ${i#truebed_} ; done ; cd .. ; done 
#have to sort before bedtools merge
for j in DEL DUP INV INS ; do cd $j ; for i in *wIDs.bed ; do sort -k1,1 -k2,2n ${i} > sorted_${i} ; done ; cd .. ; done
#merging and printing the IDs
module load bedtools2; for j in DEL DUP INV INS ; do cd $j ; for i in sorted*.bed ; do bedtools merge -c 4 -o distinct -i ${i} > merged_${i} ; done ; cd .. ; done 
#to get all the regions of interest into 1 bed file
cat extra_genes.bed promo_regions.bed gene_regions.bed | cut -f 1-4 > allregions.bed 
#doing the intersection
for j in DEL DUP INV INS ; do cd $j ; for i in merged_sorted_* ; do bedtools intersect -a allregions.bed -b ${i} -wa -wb > intersect_allregions_${i} ; done ; cd .. ; done 

#to make track files for intersected regions
for i in intersect_allregions_merged_sorted_*INS*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"0,0,255"}' ${i} > track_${i} ; done
for i in intersect_allregions_merged_sorted_*DEL*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"255,0,0"}' ${i} > track_${i} ; done
for i in intersect_allregions_merged_sorted_*DUP*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"0,255,0"}' ${i} > track_${i} ; done
for i in intersect_allregions_merged_sorted_*INV*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"255,255,0"}' ${i} > track_${i} ; done

for j in DEL DUP INV INS ; do for i in B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8 ; do cat ${j}/track_intersect_allregions_merged_sorted_${j}_${i}_combined_wIDs.bed >> track_intersect_allregions_${i}.bed ; done ; done

```
The resulting `track_intersect_allreagions_merged_sorted*` files are used in IGV to manually score which SVs occur within the candidate regions of each NAM founder line. This information is kept in a 1/0 format where a genome either has (1) or doesn't have (0) a given indel. The raw data is kept in `Indel_Haplotypes.csv`. 

## Testing for significant association with flowering time
Data: flowering time is from the BLUP scores in the Buckler et al 2009 Days to Anthesis Supplemental Table. Indel haplotypes come from above (`Indel_Haplotypes.csv`). 

Run `FT_Indel_Analysis.Rmd`

## Expression Differences for Significant Indels
### Create a file with the TPMs of all candidate genes for each NAM founder
#### Getting gene ids from NAM orthologous to the gene ID's in B73v5

##### find NAM coordinates for B73 genes

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

##### find associated NAM gene ID for those coordinates
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

##### create a matrix for gene IDs
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
##### Make tissue specific TPM files

Need bam files of the Tassel, Anther, Ear, and V11 leaf tissues
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

##### Filter the TPM files for NAM specific gene IDs
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
Then concatenate those long format files into `allTPM.txt` _need to rm headers so not repeated_
```
head -n 1 formatted_B73_candTPM.txt > allTPM.txt
for i in formatted_*; do tail -n +2 $i >> allTPM.txt ; done
```

Add columns for Genome, Tissue, and Replicate by dividing the 5th field, `Tissue_ID`
```
echo Genome > genome.txt
echo Tissue > tissue.txt
echo Replicate > replicate.txt
tail -n +2 allTPM.txt | cut -f 5 | cut -f 1 -d _ >> genome.txt
tail -n +2 allTPM.txt | cut -f 5 | cut -f 2,3 -d _ >> tissue.txt
tail -n +2 allTPM.txt | cut -f 5 | cut -f 4 -d _ >> replicate.txt

paste allTPM.txt genome.txt tissue.txt replicate.txt >> allTPM_splitdescriptions.txt
```
### Expression analysis
Run `expression_analysis.R` to get figures and results using the input files created from the above. 

## Distance to nearest QTL marker
Downloaded the legacy snp data from https://cbsusrv04.tc.cornell.edu/users/panzea/download.aspx?filegroupid=8
```
#getting the fasta sequences out of `context_full.fasta` for each QTL marker
while read -r line; do grep -A1 -w "${line}" context_full.fasta >> d2s_QTL.fasta ; done < days2silk_QTLmarkerIDs.txt
## to see which marker didn't have a position
while read -r line; do echo ${line} ; grep -c ${line} SNPpos_agpv3.bed ; done < days2silk_QTLmarkerIDs.txt	
##to get snp positions +/- 50 bp on v3 coordinates
awk -v s=50 -v OFS="\t" '{print $1,$2-s,$3+s,$4,$6}' asi_QTL.bed
##to get the fasta from the v3 reference using the bed tools
module load bedtools
bedtools getfasta -name -fi B73_RefGen_v3.fa -bed asi_QTLv3.bed > asi_QTLv3.fasta
bedtools getfasta -name -fi B73_RefGen_v3.fa -bed d2s_QTLv3.bed > d2s_QTLv3.fasta
```
## Overlap of Significantly Associated Indels with TEs
data: `NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa`
### Get the sequence of the significant deletions
* pull out the coordinates (B73v5) for the significant deletions
 `for i in track_intersect_allregions_*; do grep -f sig_DELs $i | cut -f 1-4 |sort | uniq >>sig_DELs_coords.bed; done`
 `sort -k1,1 sig_DELs_coords.bed | uniq > uniq_sig_DELs.bed`
 _since the merging of indels left some coords exactly the same while differing labels, I manually removed duplicate coordinates_
* Use the coordinates to pull out sequence from the B73v5 assembly
```
module load bedtools2
bedtools getfasta -name -fi ref_genomes/B73.PLATINUM.pseudomolecules-v1.fasta -bed uniq_sig_DELs.bed > sig_DELs.fasta
```
### Blast sequences against TE library
`nova:/work/LAS/mhufford-lab/B73/NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa`
```
module load blast-plus/2.7.1-py2-vvbzyor
blastn -max_hsps 1 -subject sig_DELs.fasta -query NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa -out out_sig_DELs_TELibrary_withlengths -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs qlen slen"
bash scripts/filterblastn.sh out_sig_DELs_TELibrary_withlengths
awk -v OFS='\t' '{print $1,$2,$6/$10, $6/$11}' out_filtered_out_sig_DELs_TELibrary_withlengths | sort -k 2,2> besthits_percov_sigDELs_TEs.txt
#manually add in a header
```
