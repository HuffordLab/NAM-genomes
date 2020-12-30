# Tandem repeat identification

1. scripts for detecting [telomere](tandemrepeat-quantification/scripts/telomere.sh):
   
   Split the pseudomolecule fasta file (using the split-fasta python package) for each NAM line to individual chromosome sequences in FASTA format : 
  
   splitfasta pseudomolecule.fasta

   Telomere script: This script takes as input the chromosome sequences in FASTA format and generates output files in txt format containing telomeric coordinates      for the small and long end of the input chromosome sequence. The script should be run in a directory containing the input FASTA files. 
   
   sh telomere.sh

2. scripts for sub-telomere ([1](tandemrepeat-quantification/scripts/subtelomeres_clusters.sh), [2](tandemrepeat-quantification/scripts/subtelomeres_blast.sh)) detection.
3. [TE quantification](tandemrepeat-quantification/scripts/TE_quantification.sh) and [repeat estimation](tandemrepeat-quantification/scripts/illumina_repeat_estimation.sh) scripts
4. [Plotting](tandemrepeat-quantification/scripts/stacked_barplot.R)
