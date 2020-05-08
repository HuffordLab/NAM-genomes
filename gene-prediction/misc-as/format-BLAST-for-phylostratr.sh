#!/bin/bash
blastin=$1
taxaid=$(echo ${blastin%.*} | rev | cut -f 1 -d "_" | rev)
blastout=${blastin%.*}.tab
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore\tstaxid" > $blastout
awk -v x=${taxaid} 'BEGIN{OFS=FS="\t"}{print $1,$2,$7,$8,$9, $10, $11, $13, x}' ${blastin} >> $blastout

