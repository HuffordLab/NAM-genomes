# GATK SNP calling

_**all scripts are located in [`snp-calling`](structural-variation/snp-calling) directory.**_

1. Input files were downloaded from CyVerse using `getFiles.sh` script. The files that were large (ran without barcodes), were split to smaller chunks for quick processing using `split-fastq.sh` script.
2. The genome file was downloaded from CyVerse (B73.v5) and processed using `gatk-prepare-reference.sh` to create all the necessary files for running GATK pipeline.
3. The fastq files were mapped to B73.v5 and processed using `process-fastq.sh` script. Briefly, this script:
   * converts unmapped fastq to bam `FastqToSam`
   * runs Picard `MarkIlluminaAdapters`
   * converts bam back to fastq `SamToFastq`
   * Maps fastq files to B73.v5 using `bwa mem` and converts to bam file using `samtools`
   * merged unmapped reads with mapped reads using `MergeBamAlignment`
   * runs picard's `MarkDuplicates` to mark optical duplicates.
4. As a final step of processing, using `run-add-readgroups.sh` correct read groups were added to the bam files and indexed.
5. GATK was run on 1Mb intervals, using the script `gatkcmds-round-1.sh` and the intervals file `B73.PLATINUM.pseudomolecules-v1_1mb_coords.bed`, the commands were generated and was run on the cluster creating slurm job submission script using GNU parallel.
6. Once the VCF files were generated (2,813 total), they were gathered and processed to filter and retain very high quality SNPs only, using the script `gatk-process.sh`
7. The bam files were recalibrated using the filtered first round SNP files using `gatk-bsqr.sh`
8. Using the recalibrated BAM files, GATK was ran again on 1Mb intervals, using the script `gatkcmds-round-2.sh` and the intervals file `B73.PLATINUM.pseudomolecules-v1_1mb_coords.bed`, the commands were generated and was run on the cluster creating slurm job submission script using GNU parallel.
9. The final files were filtered again using the `gatk-process.sh` script again.
10. The final files were uploaded to CyVerse
