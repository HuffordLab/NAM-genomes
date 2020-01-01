# Gene Prediction

## 1. Evidence Based Gene Prediction

Overview of steps:

The steps for evidence based prediction are as follows:

1. Map RNAseq reads to the genome
2. Assemble transcripts using various transcript assemblers
3. Run Mikado to consolidate transcripts and run programs to identify splice junctions, ORFs and full length matches to SwissProt plants proteins.
4. Pick transcripts for each locus and finalize evidence based annotations.
5. Post-processing of evidence based predictions.

Steps in detail:

### 1. Map RNAseq reads:

Reads were mapped using STAR mapping program as follows:

```bash
# index
runSTARmap_index.sh NAM.pseudomolecules.fasta
# first round mapping
for read1 in *R1.fq.gz; do
  runSTARmap_round1.sh NAM.pseudomolecules.fasta ${read1}
done
# consolidate splice info
awk -f sjCollapseSamples.awk *_SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > SJ.all
# second round mapping
for read1 in *R1.fq.gz; do
  runSTARmap_round2.sh NAM.pseudomolecules.fasta ${read1}
done
# merge BAM files
merge-bam-files.sh NAM
```

### 2. Transcript Assembly

Five different transcript assemblers were ran on the merged BAM file (`merged_NAM.bam`) as follows:


```bash
runStrawberry.sh merged_NAM.bam
runStringtie.sh merged_NAM.bam
runTrinity.sh merged_NAM.bam
runClass2.sh merged_NAM.bam
runCufflinks.sh merged_NAM.bam
```

Since Trinity generates a fasta file, we will use GMAP to map them and create a GFF3 file.

```bash
runGMAPcdna.sh NAM.pseudomolecules.fasta
```

We will also need splice junctions, that we can generate using PortCullis:

```bash
runPortcullis.sh merged_NAM.bam
```

### 3. Consolidate transcripts and prepare files

Here, we will use all the transcript assemblies generated in the previous step, along with splice junctions and pool them together. The redundant transcripts are also removed. Using consolidated transcripts, we will find full length trasncritps by BLAST searching against SwissProt proteins, and detect ORFs using TransDecoder.

This is all accomplished using this [script](copy-files-for-mikado.sh):

```bash
copy-files-for-mikado.sh NAM
```

### 4. Finalize evidence-based predictions

As a last step, we will run Mikado to finalize the predictions. It is done using this [script](finalize-files-and-pick.sh):


```bash
finalize-files-and-pick.sh NAM
```

### 5. Post processing of evidence-based predictions

Kapeel's section
