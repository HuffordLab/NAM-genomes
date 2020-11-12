test_bed_1=$1
test_bed_2=$2
genome=$3
output_file=$4

intersect_size=$(bedtools intersect -a $test_bed_1 -b $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
test_bed_2_size=$(bedtools merge -i $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
test_bed_1_size=$(bedtools merge -i $test_bed_1 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
genome_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $genome)
((test_bed_1_subtract_test_bed_2=$test_bed_1_size-$intersect_size))
((test_bed_2_subtract_test_bed_1=$test_bed_2_size-$intersect_size))
((genome_subtract_test_beds=$genome_size-$intersect_size-$test_bed_1_subtract_test_bed_2-$test_bed_2_subtract_test_bed_1))
echo -e "Intersect\t$intersect_size" > $output_file
echo -e "Test_bed_1_subtract_test_bed_2\t$test_bed_1_subtract_test_bed_2" >> $output_file
echo -e "Test_bed_2_subtract_test_bed_1\t$test_bed_2_subtract_test_bed_1" >> $output_file
echo -e "Genome_subtract_test_beds\t$genome_subtract_test_beds" >> $output_file
