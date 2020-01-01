#!/bin/bash
set -x
# generate List file
nam=$1
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Stringtie/merged_${nam}_stringtie.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_Stringtie.gtf
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Strawberry/strawberry-merged_${nam}.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_Strawberry.gtf
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Cufflinks/merged_${nam}/transcripts.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_Cufflinks.gtf
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Class2/merged_${nam}_class.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_class2.gtf
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/Trinity/trinity_out_dir/Trinity-GG.fasta /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/${nam}_TrinityGG.fasta
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado
echo -e "$(realpath ${nam}_Stringtie.gtf)\tSt\tTrue" > list.txt
echo -e "$(realpath ${nam}_Strawberry.gtf)\tSw\tTrue" >> list.txt
echo -e "$(realpath ${nam}_Cufflinks.gtf)\tCf\tTrue" >> list.txt
echo -e "$(realpath ${nam}_class2.gtf)\tCl\tTrue" >> list.txt
echo -e "$(realpath ${nam}_TrinityGG-mapped.gff3)\tTr\tFalse" >> list.txt
i=0;
list=$(realpath list.txt)
# get files required
junctions=$(find $(pwd) -name "portcullis_filtered.pass.junctions.bed")
genome=$(find $(pwd) -name "*pseudomolecules-v*.fasta")
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado
# run configure
mikado configure \
   --list $list \
   --reference $genome \
   --mode permissive \
   --scoring plants.yaml \
   --copy-scoring plants.yaml \
   --junctions $junctions configuration.yaml
# edit config files
cp /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/plants.yaml ./
sed -i 's/procs: 1/procs: 36/g' configuration.yaml
sed -i 's/plants.yaml/\/work\/LAS\/mhufford-lab\/arnstrm\/newNAM\/analyses\/g-transcript-assembly\/'${nam}'\/mikado\/plants.yaml/g' configuration.yaml
# consolidate transcripts
mikado prepare \
   --json-conf configuration.yaml
mkdir blastjobs
cd blastjobs
ln -s ../mikado_prepared.fasta
# split fasta and run BLAST
fasta-splitter.pl --n-parts 8 mikado_prepared.fasta
unlink mikado_prepared.fasta
for f in mikado_prepared.part-?.fasta; do
   echo "/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/runBLASTx.sh $f"
done > ${nam}.cmds
makeSLURMp.py 8 ${nam}.cmds
for f in /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/uniprot-sprot_viridiplantae.fasta*; do
ln -s $f;
done
for f in *.sub; do 
 sed -i 's/j 1 --joblog/j 8 --joblog/g' $f;
 sbatch $f;
done
# run Transcdecoder
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado
echo "/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/runTransDecoder.sh" > td.cmds
makeSLURMs.py 1 td.cmds
sbatch td_0.sub

