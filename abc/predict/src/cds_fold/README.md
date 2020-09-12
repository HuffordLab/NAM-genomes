# cds_fold
Outputs the 0-fold and 4-fold sites from a gff and reference genome.

NOTE: Assumes there is one mutation per codon.

Requires `biopython`

`biopython` can be installed with:

`pip install biopython`

Usage:

`python cds_fold.py -o sitefold.txt <reference.fasta> <reference.gff> > reference.log`


-o with a valid file name is required. Information about failed sites is written to stdout.

```
usage: cds_fold.py [-h] -o OUTFILE reference gff

Generate list of 0-fold and 4-fold sites from a fasta formatted reference
genome and gff annotation file.

positional arguments:
  reference             Fasta formatted reference file
  gff                   version 3 gff file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Name of the output file
```


Output is tab delimited with three columns: the chromosome name, the one-indexed position, and a 0 or 4 indicating if site is 0-fold or 4-fold, respectively. For example:

```
#chrom pos fold
chr1 41529 0
chr1 41531 0
chr1 41533 4
chr1 41534 0
```
