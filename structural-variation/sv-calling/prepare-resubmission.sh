#!/bin/bash
nam=$1
cd $nam
mkdir -p done-samfiles;
for f in $(cat *.log | awk '$7==0' | cut -f 9 |awk '{print $2}' |sed 's/error_corrected_reads\//B73.PLATINUM.pseudomolecules-v1_/g' |sed 's/.fasta.gz/.fasta.sam/g'); do
mv $f done-samfiles/;
done
mkdir -p first-round
mv ${nam}* ./first-round/
cd first-round
cat *.log >> ../all.log
cd ..
filter_parllel_log.sh all.log first-round/*.cmds > ${nam}b.cmds
wc -l ${nam}b.cmds
makeSLURMp.py 130 ${nam}b.cmds
for f in *.sub; do sed -i 's/parallel -j 1 --joblog/parallel -j 9 --joblog/g' $f; done

