# Centromere positions of each NAM line were projected to B73 by mapping both CENH3 ChIP-seq data and genomic input data to the B73 genome. 

## 1. Read alignment
 Both CENH3 ChIP-seq data and genomic input data to the B73 genome with bwa-mem (v0.7.17).
 ```
 genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/B73.PLATINUM.pseudomolecules-v1.fasta
for file in /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.subsample.fq.gz
do
  prefix=$(basename $file | cut -f1 -d ".")
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/SE_alignment.sh $file $genome
done
```
## 2. Enrichment calculation  
ChIP enrichment was calculated by normalizing RPKM values from the ChIP data against the genomic input in 5 kbp windows with deeptools (v3.3). 
 ```
 genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/B73.PLATINUM.pseudomolecules-v1.fasta
for file in /scratch/jl03308/NAM_pancentromere/analysis/peak_call/*/*.ChIP.q20.sorted.bam
do
  prefix=$(basename $file | cut -f1 -d ".")
  input=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/${prefix}/${prefix}.q20.sorted.bam
  sh /home/jl03308/git/NAM_pancentromere/NAM_centromere_toB73/SE_alignment.sh $file $input
done
```
## 3. Island identification
Enriched islands with a ratio above 2.5 were identified and merged with a distance interval of 1 Mbp using bedtools (v2.29).
```
for file in /scratch/jl03308/NAM_pancentromere/analysis/peak_call/*/*.ChIP_Input.RPKM.bedgraph
do
  prefix=$(basename $file |cut -f1 -d ".")
  cat $file | awk '{if($4>2.5){print$0}}'|bedtools merge -i - -d 1000000 -o mean -c 4 |grep "chr"  | \
  awk -v var1=$prefix '{print var1"\t"$0}' >> /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed
done
```

Centromeres that were not mappable by CENH3 ChIP were defined as the midpoint of the largest CentC array in B73.
```
lines=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed |cut -f1 |sort |uniq)
for line in $lines
do
  for chr in chr{1..10}
  do
    #echo -e $line"\t"$chr
    loc=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed | \
    awk -v var1=$line -v var2=$chr '{if($1==var1&&$2==var2){print$0}}') 
    if [[ ! "$loc" ]]
    then
    b73_loc_start=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/B73.majorcentc.bed | grep "$chr" | \
    cut -f3,4 | awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; printf "%.0f\n", sum}')
    b73_loc_end=$(echo $b73_loc_start + 100 |bc)
    echo -e $line"\t"$chr"\t"$b73_loc_start"\t"$b73_loc_end >> /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.all.bed
    else
    loc_start=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed | \
    awk -v var1=$line -v var2=$chr '{if($1==var1&&$2==var2){print$3}}')
    loc_end=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed | \
    awk -v var1=$line -v var2=$chr '{if($1==var1&&$2==var2){print$4}}')
    echo -e $line"\t"$chr"\t"$loc_start"\t"$loc_end >> /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.all.bed
    fi
  done
done
```
The final centromere coordinates were determined by visual inspection of the ChIP-seq peaks in IGV (v2.8).
