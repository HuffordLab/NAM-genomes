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
# GL15 as an Identified Candidate
## Getting the NAM Coordinates for GL15
GL15 V3 ID is actually: `GRMZM2G160730` and V5 ID is `Zm00001eb387280`
pan-gene id is `Pan_gene_19318`
```
grep GRMZM2G160730 Li2016_candidates/B73v3_B73v5_liftover_genemodel_CDS_xref_shortened.txt
#in R from the Li2016_Candidate_Analysis.R
> grep("Zm00001eb387280", geneIDkey$B73_AltID2)
[1] 92361
geneIDkey[92361,] %>% t() 
geneIDkey[92361,] %>% write_csv("../../gl15_nam_IDs.csv") 
tail -n 1 gl15_nam_IDs.csv | tr ',' '\n' | sort | uniq >> grep_ID.txt
```
On Nova to get the actual coordinates for each NAM
```
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/t-pangene-matrix/counts-matrix/counts-padded
for i in *.txt ; do awk -v OFS='\t' '$1 ~/Pan_gene_19318/ {print $0}' $i ; done 
#The above is saved in GL15NAMcoordinates.txt

awk -v OFS='\t' '{print $5,$6,$7,$3"_"$1"_"$2,$8}' GL15NAMcoordinates.txt > GL15NAMcoordinates.bed
```
Use nova coordinates to get the fasta for GL15 by each NAM
Need to make each genome's coordinates a different file
```
while read -r line ; do 
	n=$(echo $line | cut -f 4 | cut -f 3 -d _)
	echo $line > ${n}_GL15.bed 
done < GL15NAMcoordinates.bed
for i in *GL15.bed ; do sed -i "s/ /\t/g" $i ; done
module load bedtools2
for i in GL15_additional_work/*GL15.bed ; do 
	n=$(echo ${i#GL15_additional_work/} | cut -f 1 -d _ )
	bedtools getfasta -name -fi ref_genomes/${n}.pseudomolecules*.fasta -bed $i > GL15_additional_work/${n}_GL15.fasta 
done
cat *GL15.fasta > GL15_NAMseq.fasta 
```
Run muscle alignment on all the GL15 fastas
```
#muscle/3.8.1551 is the module used
module load muscle
muscle -in GL15_NAMseq.fasta -fastaout GL15_alignment.fasta -htmlout GL15_alignment.html -log GL15_NAMseq_aln.log
```
* Look for the insertion sequence *
Find the coordinates for the insertion sequence on P39, blast against TE library
Looking at the B73 coordinates for the insertions and gene,
I estimate the insertions of interest are between 3000 and 500 bp upstream of the gene coordinates
```
for i in *_GL15.bed ; do
awk -v OFS='\t' '{print $1, $2-3000, $2-500, $4, $5}' $i > ${i%.bed}_promoter.bed 
done
for i in GL15_additional_work/*GL15_promoter.bed ; do 
	n=$(echo ${i#GL15_additional_work/} | cut -f 1 -d _ )
	bedtools getfasta -name -fi ref_genomes/${n}.pseudomolecules*.fasta -bed $i > GL15_additional_work/${n}_GL15_promoter.fasta 
done
cat *GL15_promoter.fasta > GL15_promoter_NAMseq.fasta 
muscle -in GL15_promoter_NAMseq.fasta -fastaout GL15_promoter_alignment.fasta -htmlout GL15_promoter_alignment.html -log GL15_promoter_NAMseq_aln.log
```
*Maggie is using COGE to find the exact coordinates of the insertion* 
```
#On CML228's coordinate system: chr9 105190700-105191100
bedtools getfasta -name -fi ref_genomes/CML228.pseudomolecules*.fasta -bed GL15ins_coords_CML228.bed > CML228_GL15ins.fasta
module load blast-plus/2.7.1-py2-vvbzyor
blastn -subject CML228_GL15ins.fasta -query NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa -out out_CML228GL15ins_withlengths -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs qlen slen"
blastn -subject CML228_GL15ins.fasta -query ref_genomes/B73.PLATINUM.pseudomolecules-v1.fasta -out out_CML228GL15ins_B73assembly.txt -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs qlen slen"
```

# Li et al 2016 Candidate Analysis
## 1. Get data set up
- convert Li Candidate names to V5, filter to only those found in US-NAM (+any combination)
`Li2016_suppdata8_filtered.txt` and `B73v3_B73v5_liftover_genemodel_CDS_xref_shortened.txt`
```
grep -f Li2016-v3ids.txt B73v3_B73v5_liftover_genemodel_CDS_xref_shortened.txt >Li2016-v5ids.txt
```
- Get GWAS results
`https://iastate.box.com/s/kbmdpkydae0lvq789wmwph39l0hybi47`
```
T2	Days_To_Silk
T3	ASI
T30	Days_To_Anthesis
```

## 2. GWAS, candidate overlap and permutation
Look for the overlap of Li candidates and GWAS hits from NAM GWAS. 
Use the same window as before, gene + 5kb promoter region
Permutation test: do we see more GWAS hit overlap with these candidates compared to whole genome

Got the coordinates from grepping the V3 ID's of candidate lists against the `gene5KB_coords.bed` file
Use `intersect.sh` on condo plus significant SNP bed files from the R markdown
Also need `input.fofn` which is a File Of File Names to run as an array job submission
```
grep "^coordinates_" *.out | datamash -s crosstab 1,2 sum 5  > col5 #median SNP hits per gene
grep "^coordinates_" *.out | datamash -s crosstab 1,2 sum 4 > col4 #number of genes with at least 1 SNP hit
grep "^coordinates_" *.out | datamash -s crosstab 1,2 sum 3  > col3 #total genes sampled
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3 FS $4;next}{ print $0, a[$1]}' col4 col3 > temp1 #join cross tabbed files
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3 FS $4;next}{ print $0, a[$1]}' col5 temp1 > permutation_stats.txt #final join
```

## 3. Coefficient of Variation of candidate genes, permutation test
- Calculate the coefficient of variation for each gene's expression across all lines
- How variable is expression of genes across all lines? 
- Do we see greater variation in expression of Li candidate genes vs random sample of genes?

Coefficient of Variation is calculated as the std. dev. / mean 
`https://www.statisticshowto.com/probability-and-statistics/how-to-find-a-coefficient-of-variation/`
We need the mean and standard deviation of a gene's expression across tissue and genomes. 
Current TPM files have the following headers
```
 head -n 1 B97_joined_cnts_wLengths_tpm.txt 
Length	samples	cnts_B97_R1_anther_MN03081.txt	cnts_B97_R1_anther_MN03082.txt	cnts_B97_V11_base_MN03031.txt	cnts_B97_V11_base_MN03032.txt	cnts_B97_V11_middle_MN03041.txt	cnts_B97_V11_middle_MN03042.txt	cnts_B97_V11_tip_MN03051.txt	cnts_B97_V11_tip_MN03052.txt	cnts_B97_V11_tip_MN03053.txt	cnts_B97_V18_ear_MN03071.txt	cnts_B97_V18_ear_MN03072.txcnts_B97_V18_tassel_MN03061.txt	cnts_B97_V18_tassel_MN03062.txt
```
So for a given row (gene):
grab the TPM values across the columns (tissues and reps for a single genome) for each tpm file (all genomes)
find the mean and standard deviation of those values
and report gene + mean + standard deviation (could use awk)

```
awk '{ A=0 ; V=0; for(N=1; N<=NF ; N++) A+=$N ; A/=NF ; 
	for (N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V), A/=NF }' file
```
Test file
```
 col1  col2  col3  col4  col5
    1     2     3     4     5
    2     4     6     8    10
    3     6     9    12    15
```

### Is there a relationship between the significance of the overlap SNPs and amount of variation?
Then if there is a strong trend, that would be good follow up to structural variants. 

## 4. Differential Expression between tropical and temperate lines
- Group lines by tropical and temperate lines (remove mixed lines)
- Is there a significant difference in expression between these line groups for candidate genes?
- If there is, are we more likely to find SVs around these candidates than genes in general? (Permutation test)

### Looking for SVs nearby Dof21 which was the only thing that popped in #4
Using the full genome tracks created up in the Arun's method section, `SV_v8`, IGV
`GRMZM2G162749` == `Dof21` == `Zm00001eb005380_T001`
 From the 5KB+gene window bed file `gene5KB_coords.bed`
 ```
 1	14950395	14957579	GRMZM2G162749
 ```
Looked at log2 fold change to find other candidates: ZCN10, ZCN8, ZMCCT10

