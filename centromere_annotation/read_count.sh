for file in /scratch/jl03308/NAM_pancentromere/rawdata/ChIP_new/*.fq.gz_trimming_report.txt
do
  prefix=$(basename $file |cut -f1 -d ".")
  val=$(cat $file |grep "Total reads processed:" |cut -f2 -d ":"| sed "s/^[ \t]*//")
  echo -e $prefix"\t"$val
done
