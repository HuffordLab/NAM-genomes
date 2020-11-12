#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 3GB of memory for each process
#$ -pe threads 30 -l m_mem_free=0.5G
# make sure that there is at least 20G of tmp space on the node
#$ -l tmp_free=300G

date
lines=( B73 CML247 CML333 HP301 Ki3 M37W NC350 Oh7b Tzi8 B73_AB10 CML103 CML277 CML52 IL14H Ky21 Mo18W NC358 P39 B97 CML228 CML322 CML69 Ki11 M162W MS71 Oh43 Tx303 )

for i in "${lines[@]}"
do
SPECIE=$i
echo "$SPECIE"

python ../TEsorter.py ${nam_specie}.merged.cds.fasta -eval 1e-6 -p 30

grep -v "^#" ${SPECIE}.maker.cds.fasta.rexdb.cls.tsv | cut -f1 | sort | uniq | cut -f1 -d "_" | sort | uniq > TE-genes.txt

grep -Fvf TE-genes.txt ${nam_specie}.merged.gff | awk '$3 !~ /superlocus/' > temp.gff3

grep -Ff TE-genes.txt ${nam_specie}.merged.gff | awk '$3 !~ /superlocus/' > temp.other.gff3

gt gff3 -sort -tidy -retainids temp.gff3  > ${SPECIE}.TErmv.gff3
gt gff3 -sort -tidy -retainids temp.other.gff3  > ${SPECIE}.withTE.gff3

mkdir $SPECIE

mv *gff3 *rexdb* TE-genes.txt $SPECIE
done
date                                                                                                                           
                                                                                                            1,1           Top
