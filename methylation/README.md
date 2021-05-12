# Methods for methylation and accessible chromatin region analyses:


## Steps

1. Download the genome fasta file and annotation files and creates the directory structure described in the [script](scripts/download_genomes.sh), using the [file](assets/NAM_founder_genome_list.txt).

2. Process the NAM genomes

    * alters the chromosome names in the NAM genome fasta files
    * creates a bedtools genome file for use in bedtools
    * creates bed files of 100 bp windows, sliding 20 bp, genome-wide
    * identifies all contiguous sequences (contigs) across the chromosomes
    * this script works in conjunction with the directory structure created by download_genomes.sh

  Scripts:

  - [submit_process_NAM_genomes.sh](scripts/submit_process_NAM_genomes.sh)
  - [process_NAM_genomes.sh](scripts/process_NAM_genomes.sh)
  - [findgaps.py](scripts/findgaps.py)

  Additional files:

  - [NAM_founder_genome_list.txt](NAM_founder_genome_list.txt)

  Additional software:

  - bioawk 1.0
  - Biopython 1.68
  - BEDTools 2.28.0

3. Process the NAM gene annotations

    * converts gff gene annotation to bed file
    * extracts gene information
    * change chromosome names
    * identify the chromosomes and contigs that have annotated genes and those that do not have annotated genes

  Scripts:

  - process_gene_annotations.sh (provided here on github)

  Additional files:

  - NAM_founder_genome_list.txt (provided here on github)

  Additional software:

  - BEDTools 2.28.0

4. Create Bowtie2 index files for all the NAM founders

    * creates a directory for each NAM founder which contains the index files for aligning reads with Bowtie2
    * Creates a separate script for each NAM founder and submits it separately

  Scripts:

  - create_bowtie2_indices.NAM.v1.sh (provided here on github)

  Additional files:

  - NAM_founder_genome_list.txt (provided here on github)

  Additional software:

  - Bowtie2 2.3.5.1

5. Create files of all cytosine coordinates in the NAM genomes

    * A separate script is created for each NAM founder and submitted separately
    * Split chromosomes into separate files for running fuznuc
    * Extract CG, CHG, and CHH coordinates from the genome fasta file
    * concatenate the split chromosomes back together to create single genome-wide file for each sequence context and NAM founder
    * this concatenated file is sorted and ready to be joined to the methylome data (which is done in a separate script)

  Scripts:

  - Extract_all_Cs_NAM_genomes.v4.sh (provided here on github)

  Additional files:

  - NAM_founder_genome_list.txt

  Additional software:

  - pyfaidx 0.5.5.1
  - EMBOSS 6.6.0

6. Processing the methylome CGmap files

   * A separate script is created for each tissue and submitted separately.
   * the CGmap files are listed in field 5 of the NAM_methylome_list_merge.v2.txt file.
   * These CGmap files are located in the CGmap_home directory
   * Convert CGmap files to allc bed files
   * correct scaffold names to match the reference genome
   * Merge the allc file biological replicates into single files
   * convert the allc file to a bed file
   * split the bed file by sequence contexts
   * Create the files containing coverage for all Cs in the genome, split into the three contexts, and also CHGCGN
   * Create coverage files for the three contexts

  Scripts:

  - convert_CGmap_to_allc_then_merge_then_convert_bed.v2.sh (provided here on github)

  Additional files:

  - NAM_methylome_list_merge.v2.txt (provided here on github)
  - CGmap files (generated from earlier steps in the methods)

  Additional software:

  - methylpy 1.3.2
  - BEDTools 2.28.0
  - kenttools 1.04.00

7. Call the unmethylated regions (UMRs)
    * Each NAM founder, mapped to both B73 and to its cognate reference genome, are submitted as separate scripts
    * submit_UMR_script.v2.sh is used for submitting the separate scripts
    * three variables are assigned in submit_UMR_script.v2.sh: $GENOME, $out_dir, and $prefix
    * the UMR_algorithm_all_contexts_v1.3.NAM_founders.sh scripts calls the UMRs in each sample
    * UMRs are called in both the CHG and CHG/CGN contexts, then merged together.
    * UMRs are split into those smaller than and larger than 150 bp
    * output file: ${out_dir}/${prefix}.final_UMRs.chg.chgcgn.genic-prox-dist.${size}.bed
    * output file description:
```
      field1: chrom
      field2: UMR start
      field3: UMR end
      field4: methylation count
      field5: read coverage count
      field6: informative cytosine count
      field7: percent methylation
      field8: distance of UMR to nearest annotated gene
      field9: "genic", "proximal", or "distal" designation
```

  Scripts:

  - submit_UMR_script.v2.sh (provided here on github)
  - UMR_algorithm_all_contexts_v1.3.NAM_founders.sh (provided here on github)

  Additional files:

  - NAM_methylome_list_merge.v2.txt (provided here on github)
  - allc.bed files generated in step 6
  - files of total cytosine coverage, generated in step 6
  - 100 bp sliding window bed files generated in step 2
  - bedtools genome files generated in step 2
  - list of contigs lacking genes, generated in step 2
  - bed files of genome N-gaps, generated in step 2
  - bed files of gene annotations, generated in step 3

  Additional software:

  - BEDTools 2.28.0

8. Clean up the UMR files
    * Remove the blacklisted regions from the previously called UMRs
    * Rename the UMRs to the standardized naming scheme
    * Convert the files to the formal BED format

  Scripts:

  - remove_blacklist.change_filenames.convert_bed.NAM_UMRs.v1.sh (provided here on github)

  Additional files:

  - methylomes_table_cross_ref.txt
  - UMR files generated in step 7
  - Blacklist files

  Additional software:

  - BEDTools 2.28.0


9. Calling differentially methylated regions (DMRs) and conserved UMRs

    * A separate script is generated and submitted for each NAM line
    * The coordinates of UMRs in B73 are split into quarters
    * The mCHG data from allc files in each NAM line is mapped to the B73 UMR quarters
    * Quarters meeting the mapping criteria are kept
    *  B73 UMRs with all four quartiles below 20% mCHG or above 50% mCHG are kept, and designated conserved hypomethylated, or DMR, respectively

  Scripts:

  - extract_NAM_DMRs.v1.sh (provided here on github)

  Additional files:

  - NAM_mapping_list.txt (provided here on github)
  - 7-column UMR bed file generated in step 7
  - meth_*.ref_B73.allc.bed.chg file generated in step 6

  Additional software:

  - BEDTools 2.29.2


10. The initial ATAC-seq read processing
    * Run each ATAC rep (fastq file) through fastp in paired-end mode

  Scripts:

  - ATAC_initial_read_process.NAM.v1.sh (provided here on github)

  Additional files:

  - genome_mapping_combinations.txt (provided here on github)

  Additional software:

  - fastp 0.20.0

11. Map the ATAC-seq reads with bowtie2
    * Align each ATAC-seq replicate to the B73 reference genome and also the cognate reference genome
    * Each replicate is submitted as a separate script

  Scripts:

  - ATAC_map_reads.NAM.v1.sh (provided here on github)

  Additional files:

  - genome_mapping_combinations.txt (provided here on github)

  Additional software:

  - Bowtie2 2.3.5.1

12. Filter the aligned ATAC-seq reads
    * This script creates and submits a separate script for each NAM line, mapped to B73 reference and its cognate reference
    * Convert `sam` file to `bam` file
    * remove duplicates
    * filter by `MAPQ >= 30`
    * Create 1 bp resolution bigwig file

  Scripts:

  - ATAC_process_bams_v3.sh (provided here on github)

  Additional files:

  - genome_mapping_combinations.txt (provided here on github)

  Additional software:

  - deepTools 3.3.1
  - SAMtools 1.10
  - sambamba 0.7.1

13. Call accessible chromatin regions (ACRs) with MACS2 and quantify their overlap with UMRs
    * Call peaks on ATAC data with MACS2 using paired-end mode, q-value cutoff of 0.005, and effective genome size of 1.8e+9 bp
    * remove ACRs overlapping blacklisted N-gap regions
    * designate ACR as genic, proximal, or distal, depending on proximity to nearest annotated gene
    * calculate the percent overlap of ACRs with UMRs

  Scripts:

  - macs2_call_acrs.nam_lines.v5.sh (provided here on github)

  Additional files:

  - genome_mapping_combinations.txt (provided here on github)
  - UMR bed files generated in step 8
  - ATAC-seq bam files generated in step 11
  - list of contigs lacking genes, generated in step 2
  - bed files of genome N-gaps, generated in step 2
  - bed files of gene annotations, generated in step 3

  Additional software:

  - MACS2 2.2.7.1
  - BEDTools 2.29.2
