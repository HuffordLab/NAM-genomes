#This script will cluster the blast hits (distance between two clusters is set to 15kb) to identify the start and stop boundaries of subtelomeric repeat array at each end of chromosome of a genome
module load bedtools2

for file in *bed;do # bed files for each chromosome generated as output of the subtelomere blastscript
outname_pre=$(echo $file| cut -d'_' -f 2)
outname=$(echo $outname_pre| cut -d'.' -f 1)
bedtools cluster -i $file -d 15000 | awk '{if($4==1)print }' | head  -1 | awk '{print $2}' > $outname.start_small_arm
bedtools cluster -i $file -d 15000 | awk '{if($4==1)print }' | tail -1 | awk '{print $3}' > $outname.end_small_arm
paste $outname.start_small_arm $outname.end_small_arm > $outname.small_arm_coords.txt
last_cluster=$(bedtools cluster -i $file -d 15000 | tail -1 | awk '{print $4}') 
bedtools cluster -i $file -d 15000 | awk '$4==last_c {print}' last_c="$last_cluster" | head -1 | awk '{print $2}' > $outname.start_long_arm
bedtools cluster -i $file -d 15000 | awk '$4==last_c {print}' last_c="$last_cluster" | tail -1 | awk '{print $3}' > $outname.end_long_arm
paste $outname.start_long_arm $outname.end_long_arm > $outname.long_arm_coords.txt
done


