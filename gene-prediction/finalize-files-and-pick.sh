#!/bin/bash
nam=$1
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/blastjobs
for xml in *xml.gz; do
   gunzip $xml;
done
python ~/gitdirs/common_analyses/BlastXMLmerge.py ../mikado.blast.xml *.xml
for xml in *xml; do
   gzip $xml;
done
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/
blastxml=mikado.blast.xml
targets=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/blastjobs/uniprot-sprot_viridiplantae.fasta
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
junctions=$(find $(pwd) -name "portcullis_filtered.pass.junctions.bed")
genome=$(find $(pwd) -name "*pseudomolecules-v*.fasta")
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado
#serialise
mikado serialise \
   --start-method spawn \
   --procs 28 \
   --blast_targets ${targets} \
   --json-conf configuration.yaml \
   --xml ${blastxml} \
   --orfs ${orfs}
#pick
mikado pick \
   --start-method spawn \
   --procs 28 \
   --json-conf configuration.yaml \
   --subloci_out mikado.subloci.gff3

mkdir -p /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado
cp mikado.loci.gff3 /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/${nam}_mikado.loci.gff3
cp mikado.loci.metrics.tsv /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/${nam}_mikado.loci.metrics.tsv
cp mikado.loci.scores.tsv /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/${nam}_mikado.loci.scores.tsv
module load cufflinks
gffread mikado.loci.gff3 -g ${genome} -t mRNA -x ${nam}_mikado.transcripts.fasta -y ${nam}_mikado.proteins.fasta
cp ${nam}_mikado.proteins.fasta /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/${nam}_mikado.proteins.fasta
cp ${nam}_mikado.transcripts.fasta /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/${nam}_mikado.transcripts.fasta
mv ${nam}_TrinityGG-mapped.gff3 /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
mv ${nam}_Stringtie.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
mv ${nam}_Strawberry.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
mv ${nam}_Cufflinks.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
mv ${nam}_class2.gtf /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
mv ${nam}_TrinityGG.fasta /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/${nam}/base-assemblies/
