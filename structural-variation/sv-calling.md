# Structural variation

Structural variation was called using [SNIFFLES](https://github.com/fritzsedlazeck/Sniffles). The scripts can be found [here](sv-calling).


## 1. Read alignment

Read alignment was performed using [NGMLR](https://github.com/philres/ngmlr) aligner using B73.v5 as reference and the error-corrected PacBio reads of NAMs as query. The folders were organized such that for each NAM a separate folder containing its error corrected reads were placed.

```bash
/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/a-corrected-reads-mapping
tree -L 1
.
├── B73
├── B73_Ab10
├── B97
├── CML103
├── CML228
├── CML247
├── CML277
├── CML322
├── CML333
├── CML52
├── CML69
├── HP301
├── IL14
├── Ki11
├── Ki3
├── Ky21
├── M162W
├── M37W
├── Mo18W
├── MS71
├── NC350
├── NC358
├── Oh43
├── Oh7b
├── P39
├── Tx303
└── Tzi8
# split the reads to small chunks:
for ec in */*error-corrected.fasta.bz2; do
  cd $(dirname $ec);
  split-zip.sh $(basename $ec);
  cd ..
done
# generate SLURM script and submit jobs
for nam in $(ls -d */); do
  ./makeCMDS.sh ${nam::-1};
  cd ${nam::-1};
  for sub in *.sub; do
    sbatch $sub;
  done
  cd ..;
done
# merge bam files and clean-up
for nam in $(ls -d */); do
  ./runSAMcleaner.sh ${nam::-1};
done
```
## 2. Calling variants using SNIFFLES


First, softlink all the alignments in a new directory.

```bash
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/a-corrected-reads-mapping
mkdir SNIFFLES_v3_20200527
find $(pwd) -name "NAM_merged_ec_reads_mapped_to_B73.bam" | awk -F"/" '{print "ln -s "$0,$9".bam"}' > softlink.sh
cd SNIFFLES_v3_20200527
bash ../softlink.sh
ls *.bam > bam.fofn
```

The first round of SNIFFLES was carried out as follows:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 1:00:00
#SBATCH -J sniffles
#SBATCH -o sniffles.o%j
#SBATCH -e sniffles.e%j
IFS=$'\n' read -d '' -r -a CMDS < bam.fofn
CMD=${CMDS[$SLURM_ARRAY_TASK_ID]}
./runSNIFFLES.sh ${CMD}
```
Organize:

```bash
mkdir -p r1_bedpe/logfiles r1_vcf/logfiles
mv *bedpe.log ./r1_bedpe/logfiles
mv *vcf.log ./r1_vcf/logfiles
mv *.vcf ./r1_vcf
mv *.bedpe ./r1_bedpe
```
Merge

```bash
cd r1_vcf
../merge-first-round-vcf.sh
```

The second round of SNIFFLES was then run as follows:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 1:00:00
#SBATCH -J sniffles
#SBATCH -o sniffles.o%j
#SBATCH -e sniffles.e%j
IFS=$'\n' read -d '' -r -a CMDS < bam.fofn
CMD=${CMDS[$SLURM_ARRAY_TASK_ID]}
./runSNIFFLES-gt.sh ${CMD}
```
Organize:

```bash
mkdir r2_vcf/logfiles
mv *vcf.log ./r2_vcf/logfiles
mv *.vcf ./r2_vcf
```
Merge

```bash
cd r2_vcf
../merge-second-round-vcf.sh
```

## 3. Calling variants using BioNano Solve


Remove all the scaffolds from the genome (fasta file)

```bash
./separate-chr.sh <nam-genome>.fasta
```

Convert fasta to cmap format


```
fasta2cmap.sh <NAM-genome>-chr-only.fa
```

Call SV's against B73

```
bionano-call-SV.sh <NAM-genome>-chr-only.cmap
```

convert to vcf format

```
smap2vcf.sh
```
