# Quality Control

Methods and scripts used for checking the quality of datasets at various stages of the project are documented here. Some methods are for the raw data inspection, while others are to check the source/line integrity, as well as the final data results quality.  

Quality inspection was done for:

1. _Illumina short reads_ generated for each NAM were inspected to make sure that sequencing was done on the right NAM line before proceeding with the PacBio sequencing and assembly.
The data was compared with the HapMap data to make sure the newly sequenced line clusters with the correct accession within the HapMap data.
RawSequence obtained from `BaseSpace` were already inspected for the sequencing quality, so it was skipped.

2. _PacBio data_ generated from multiple sequencing facility was tested using **SequelTools** which provides various quality metrics for inspecting and accessing the data quality.

3. _RNA-seq data_ downloaded from BaseSpace (which was already tested for quality, demultiplexed and trimmed of any adapter sequences) and were tested by first mapping against V4 B73 (to check the mapping percent), counts were generated against genic features provided by the V4 annotation, and count based clustering using DESeq2 was performed to test correct clustering of various tissues within the NAM genome. The RNAseq data was further tested to verify its source (accession and tissue) by calling variants using GATK recommended workflow.

4. _Assembly QC_ Finished assembled sequences (contigs/scaffolds/pseudomolecule) were used as input sequences in FASTA format to the GenomeQC pipeline to calculates different quality metrics for assessing the contiguity and gene space completness of the NAM genome assemblies. For assessing gene space completeness, Embryophyta odb9 (lineage specific profile) and maize (AUGUSTUS gene prediction species profile) were selected as the input parameters for calculating % BUSCO completeness.




## Illumina QC:

1. A clean data set containing HapMap2 reference maize SNPs was prepared
	- HapMap2 maize SNP reference data was downloaded from the panzea website in the condensed HapMap format
	- Monomorphic loci were removed using a custom Python [scripts](scripts/18_1_removeMMLhapMap.py)
	- The file was converted from the condensed HapMap format to the standard Hapmap format using a custom Python [script](scripts/18_1_standardizeHapMap.py)
	- The file was then filtered so as to remove all SNPs other than those associated with relevant NAM founders and [TIL11](scripts/18_1_hapMapNAMonly2.py)
2. A clean SNP data set was generated using raw Illumina reads from one or more newly sequenced NAM founders (modified from Li's scripts: [1](https://github.com/HuffordLab/Wang_et_al._Demography/blob/master/trim_mapping_MD/trim_pe.sh), [2](https://github.com/HuffordLab/Wang_et_al._Demography/tree/master/GATK_SNPcalling), [3](https://github.com/HuffordLab/Wang_et_al._Demography/blob/master/trim_mapping_MD/20150709_trim.mapping.MD.sh)
	- read mapping was done with BWA in a bash pipeline
	- SNP calling was done with GATK  producing a vcf file
3. The HapMap reference SNPs and newly generated NAM founder SNPs were merged producing a SNP file in [standard HapMap format](scripts/18_1_mergeOurDataWHapMap.py)
4. Percent Heterozygosity was calculated from merged SNP data set with the [script](scripts/18_1_countPerPolyLoci2.py)
5. SNPphylo was used to build a phylogenec tree from the merged SNP data set
	- The SNP file was randomly subsampled such that the SNPphylo input contains no more than 50,000 SNPS ,which is a requirement of SNPphylo, [see](scripts/18_1_subSampleHapmapSNPs.py)
	- SNPphylo was used to generate the phylogenetic tree which was visualized in FigTree

Scripts used for this parts are in the [scripts](./scripts) directory.

## PacBio QC:

[SequelTools](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03751-8) was run on PacBio data to collect QC metrics. SequelTools was run as described in the manual.

```bash
# for each NAM, collect the names of scraps and subread bam files
find $(pwd) -name "*subreads.bam"  > subreads.txt
find $(pwd) -name "*scraps.bam"  > scraps.txt
# Run SequeTools with default options
bash SequelTools.sh -t Q -u subreads.txt -c scraps.txt
```
The full Excel Sheet with summary stats for each SMRTcell is [available here](assets/PacBio_SequelTools_stats.xlsx).

### Summary Stats

| NAM     | Gsize (Gb) | Facility | SMRT cells | Mean SR read length | Median SR read length | Mean N50 | Depth of coverage (SR) |
|:--------|-----------:|:---------|-----------:|--------------------:|----------------------:|---------:|-----------------------:|
| CML277  | 2.337      | Arizona  | 27         | 8840.81             | 6137.0                | 15682.2  | 70.7                   |
| KY21    | 2.286      | OSU      | 21         | 11824.67            | 9602.0                | 18338.8  | 68.7                   |
| Tzi8    | 2.346      | Arizona  | 20         | 11315.10            | 8717.0                | 18732.2  | 66.4                   |
| CML247  | 2.351      | OSU      | 22         | 10721.71            | 8001.0                | 17945.6  | 68.3                   |
| IL14H   | 2.102      | UGA      | 26         | 10795.13            | 8094.0                | 18152.4  | 85.5                   |
| Oh7b    | 2.203      | BYU      | 26         | 12302.99            | 10963.0               | 18828.7  | 70.6                   |
| HP301   | 2.208      | UGA      | 22         | 12721.05            | 10033.5               | 20426.8  | 73.8                   |
| M162W   | 2.238      | OSU      | 22         | 13293.20            | 11109.5               | 20897.4  | 76.2                   |
| Ki3     | 2.348      | UGA      | 22         | 12870.39            | 11029.5               | 19749.3  | 63.5                   |
| CML228  | 2.364      | BYU      | 28         | 11296.39            | 9103.0                | 17818.4  | 68.3                   |
| NC358   | 2.272      | Arizona  | 21         | 13711.53            | 11020.0               | 21816.6  | 71.2                   |
| CML103  | 2.309      | UGA      | 23         | 11583.03            | 9071.0                | 18756.8  | 72.2                   |
| CML52   | 2.332      | BYU      | 26         | 11234.27            | 9013.0                | 17672.0  | 69.1                   |
| OH43    | 2.222      | OSU      | 17         | 13953.02            | 11157.0               | 22216.4  | 69.9                   |
| CML069  | 2.498      | BYU      | 22         | 13392.14            | 9816.0                | 22256.7  | 64.4                   |
| CML322  | 2.305      | OSU      | 18         | 13939.05            | 10760.0               | 22004.8  | 65.6                   |
| M37W    | 2.330      | UGA      | 22         | 15110.58            | 12032.0               | 24461.1  | 65.9                   |
| MO18W   | 2.328      | OSU      | 22         | 13742.10            | 10022.0               | 23843.7  | 85.2                   |
| MS71    | 2.192      | Arizona  | 22         | 15387.96            | 11909.5               | 25324.4  | 76.1                   |
| P39     | 2.091      | BYU      | 19         | 16060.80            | 11951.0               | 26507.0  | 81.3                   |
| CML333  | 2.401      | BYU      | 21         | 15666.66            | 11660.0               | 26361.6  | 63.2                   |
| B73AB10 | 2.300      | UGA      | 22         | 15965.76            | 13083.0               | 26462.7  | 62.7                   |
| Tx303   | 2.346      | Arizona  | 20         | 16180.25            | 11968.0               | 27810.9  | 71.1                   |
| B73     | 2.300      | UGA      | 18         | 17246.77            | 12437.5               | 29671.6  | 82.9                   |
| KI11    | 2.381      | UGA      | 18         | 16579.53            | 11814.5               | 29169.9  | 65.9                   |
| B97     | 2.284      | Arizona  | 23         | 16339.04            | 13496.0               | 29428.0  | 71.3                   |
| NC350   | 2.394      | UGA      | 19         | 18654.96            | 14166.0               | 31717.6  | 65.9                   |



## RNA-Seq QC:

### General Data Quality

RNA-seq data downloaded from BaseSpace were mapped to B73.v4 genome using STAR (2-pass mapping), counts were gathered using the v4 annotation in gff3 format using `featureCounts`, followed by clustering using DESeq2. Methods are as follows:

Steps:
1. soft-link all the input-files creating directory for each NAM: [`soft-link-rnaseq-files.sh`](scripts/soft-link-rnaseq-files.sh)
2. Download B73.v4 and index the genome using STAR: [`runSTARmap_index.sh`](scripts/runSTARmap_index.sh)
3. Create commands for each NAM directory: [`star-cmds-gen.sh`](scripts/star-cmds-gen.sh), this will generate commands for
     * Map raw reads to the indexed genome using STAR (1st pass): [`runSTARmap_round1.sh`](scripts/runSTARmap_round1.sh)
		 * Aggregate splice sites to generate a file for reliable splice sites found in the genome: [`sjCollapseSamples.awk`](scripts/sjCollapseSamples.awk)
		 * Map raw reads to the indexed genome again, but using the splice sites to recalibrate alignment and generate a bam file: [`runSTARmap_round2.sh`](scripts/runSTARmap_round2.sh)
		 * Run the `featureCoutns` from subread package to get counts for each gene present in the annotation: [`runSubreads.sh`](scripts/runSubreads.sh). It also runs `multiqc` for gathering quality metrics.
		 * Clean the counts file and generate R script process it with DESeq2 package using: [`cleanSubreadCounts.sh`](scripts/cleanSubreadCounts.sh). The base R-script is: [`DESeq2.R`](scripts/DESeq2.R), which needs to be run to create plots and clustering data.
		 * Gather summary stats using scripts [`tabulator.sh`](scripts/tabulator.sh) and [`getCoverageStats.sh`](scripts/getCoverageStats.sh)
4. The plots and tables were imported to Excel and PowerPoint to evaluating the quality. The files are available here: [RNA_seq_v5.2.pptx](assets/RNA_seq_v5.2.pptx) and [RNAseq_inventory_V5.2.xlsx](assets/RNAseq_inventory_V5.2.xlsx).

### Accession and Tissue integrity

Each RNAseq library was mapped to B73.v4 and SNPs were called using GATK. The SNPs were then clustered to generate a phylogenetic tree, compared to the NAM founder phylogenetic tree (from HapMap2) to detect any sample/tissue/accession mislabeling or switching. The methods/scripts are listed below:

 1. The bam file generated from the previous step was processed using GATK and Picard toolkit using the script [`process-rnaseq-bams.sh`](scripts/process-rnaseq-bams.sh). Briefly, this scripts adds read-group information, marks duplicate reads in the mapped bam file, and splits N cigar string in mapped reads.
 2. Run GATK `HaplotypeCaller` on the processed bam files. To speed it up, the commands were generated using intervals, processing 3Mb regions at a time, simultaneously, running them on the cluster. The script used: [`gatkcmds-gen.sh`](scripts/gatkcmds-gen.sh). The jobs script wrapper [`makeSLURMp.py`](https://github.com/ISUgenomics/common_scripts/blob/master/makeSLURMp.py) was used.
 3. The VCF files for each region were gathered and filtered to generate a single VCF file.
 4. SNPphylo was used to contruct phylogenetic tree and the tree was inspected to check if the tissue samples from each accession clustered correctly. Any irregularities were flagged and removed from the final analyses (samples that could be easily obtained were re-sequenced). Script [`snp-phylo-run.sh`](scripts/snp-phylo-run.sh).


## Assembly QC

```bash
sh pipeline.sh genome.fasta output_assembly_stats output_buscometrics embryophyta_odb9 maize 32
```
