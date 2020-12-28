# Genome Assembly

Raw PacBio reads (CCS) were error corrected using `FALCON` and assembled using `canu`. The settings used were as follows:

## FALCON error correction

Settings used for error-correction were as follows:

```
#### Input options
[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=streamed-median
target=pre-assembly
skip_checks=False
LA4Falcon_preload=false

#### Data Partitioning
pa_DBsplit_option= -x500 -s400
ovlp_DBsplit_option= -x500 -s400

#### Repeat Masking
pa_HPCTANmask_option=
#no-op repmask param set
pa_REPmask_code=0,300;0,300;0,300

####Pre-assembly
# adjust to your genome size
genome_size =  2400000000
seed_coverage = 45  #adjusted according to datasets
length_cutoff = -1
pa_HPCdaligner_option= -k14 -e0.75 -s100 -l3000 -h240 -w8 -H10310
pa_daligner_option= -k18 -e0.80 -l1000 -h256 -w8 -s100
falcon_sense_option= --min_idt 0.70 --min_cov 2 --max_n_read 200
falcon_sense_greedy=False
```


## CANU assembly

`canu` assembler was run on the corrected reads as follows:

```
canu -assemble \
   -p ${FILE} \
   -d ${FILE}-out \
   genomeSize=2.4g \
   -pacbio-corrected ${FILE}.corrected.fasta \
   useGrid=true \
   minReadLength=15000 \
   minOverlapLength=3000 \
   correctedErrorRate=0.065 \
   corMhapSensitivity=normal \
   ovlMerThreshold=500 \
   utgOvlMerThreshold=200 \
   cnsMemory=12 \
   cnsThreads=12 \
   merylThreads=24 \
   batThreads=24 \
   gfaThreads=24 \
```



# Hybrid Assembly using BioNano Solve

BioNano Solve was run as suggested by BioNano [documentation](https://bionanogenomics.com/wp-content/uploads/2018/04/30205-Guidelines-for-Running-Bionano-Solve-Pipeline-on-Command-Line.pdf).
