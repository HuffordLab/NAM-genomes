## Quality Insepection

Illumina reads generated for each NAM were insepcted to make sure that sequencing was done on the right NAM line before proceeding with the PacBio sequencing and assembly. The data was compred with the
HapMap data to make sure the newly sequenced line clusters with the correct accession within the HapMap data.

Methods overiview:

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
