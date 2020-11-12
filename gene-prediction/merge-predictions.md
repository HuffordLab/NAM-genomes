# Model refinements and QC

The Mikado and non-verlapping Braker predicted models were combined to have a Working Set(WS) and qaulity filtered to a High Confindent Set(HCS) using below steps:

## Step 1. Generate non-verlapping Braker set
Arun can you add the steps you used to generate this


## Step 2. Combine Mikado and non overlapping BRAKER models to a WS

Here use accessory script from MAKER-P tool to merge gff

```
gff3_merge -o ${nam_specie}.merged.gff <Mikado.gff> <non-overlapping.braker.gff>
```


## Step 3. Assignment and filtering by Transposable Element(TE).

Transposon element (TE) related genes were filtered using the [TEsorter](https://github.com/zhangrengang/TEsorter) tool which uses the REXdb (viridiplantae_v3.0 + metazoa_v3) database of TEs. Any TE-related genes from the WS was removed using the submit script [TEsorter.sh](/gene-prediction/cshl-scripts/TEsorter.sh)
```
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

```



## Step 4. Assign Annotation Edit Disctance (AED) score to WS using MAKER-P tool 
MAKER-P provides AED scores as a mersure of qaulity of annotations. THe score ranges from a scale of 0-1. With score 0 being models high agreement with the underlying evidence and 1 being least agreeable. We used MAKER-P to provide AED scores to the WS models that were TE removed and used a filter of AED < 0.75  to generate HCS. 

3.1 Before assigning the AED scores we used RepeatMasker within MAKER with curate MTEC library for repeatmasking the NAM genomes. The MAKER control files are documented here and the submit script [maker_repeatmask.sh](/gene-prediction/cshl-scripts/maker_repeatmask.sh).

3.2 To assign AED scores the HCS gff was passed as `model_gff` in the maker_opts.ctl file and the repeatmasked genome as `genome` input. Complete description of the control is prvided here and the maker submit script [maker_run.sh](/gene-prediction/cshl-scripts/maker_run.sh).

3.3 The AED assigned gff was curated and filtered for AED < 0.75 using MAKER accessory scripts as descrined below

```
gff3_merge -o ${nam_specie}.gene_only.gff <MAKER.gff>

quality_filter.pl -a 0.75 ${nam_specie}.gene_only.gff

```

## Step 5. Assignment and filtering by Phylogenetic context

Each HCS gene was assigned a phylogenetic assignment based on USEARCH 11.0.667 alignment of protein sequences to the longest isoform from Sorghum bicolor v3.1, Oryza sativa japonica IRGSP-1.0, Brachypodium distachyon v3.0 or Arabidopsis thaliana TAIR10. Genes having an isoform with an alignment to any of these annotations are considered conserved, and those lacking such an alignment are considered lineage specific genes

```
usearch11.0.667_i86linux32 -usearch_local <accession>.protein.fasta -db outgroups.udb -id 0.5 -uc <accession>.uc

```
The usearch output file was used to guide extracting the conserved and species specific genes from the fasta and gff files.
```
./make_t2g.pl <accession>.gff > <accession>.t2g
./get_species_specific_t2g.pl <accession>.t2g <accession>.uc > <accession>.t2g.ss
./conserved_filter.pl <accession>.conserved <accession>.t2g <accession>.protein.fasta <accession>.cds.fasta <accession>.gff
./non-conserved_filter.pl <accession>.specific <accession>.t2g.ss <accession>.protein.fasta <accession>.cds.fasta <accession>.gff
```

