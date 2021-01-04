
ml BEDTools/2.29.2-GCC-8.3.0
for file in /scratch/jl03308/NAM_pancentromere/analysis/peak_call/*/*.ChIP_Input.RPKM.bedgraph
do
  prefix=$(basename $file |cut -f1 -d ".")
  cat $file | awk '{if($4>2.5){print$0}}'|bedtools merge -i - -d 1000000 -o mean -c 4 |grep "chr"  | \
  awk -v var1=$prefix '{print var1"\t"$0}' >> /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM_B73.centromere.bed
done

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
