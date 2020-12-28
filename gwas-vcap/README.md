# GWAS SVs
The number of RILs is 4,027 and the number of markers is 71,196.
Mixed linear model based association analysis.

### Output file format
.mlma (columns are chromosome, Markers, physical position, reference allele (the coded effect allele),
the other allele, frequency of the reference allele, SNP effect, standard error and p-value).

The significant QTLs are defined by -log10(p-value) > FDR (0.05).
```
library("IHW")
get_bh_threshold(GWAS.Results$p[!is.na(GWAS.Results$p)], 0.05)
```

For a few of traits, FDR threshold is not there, 5 is used.

36 output files are for 36 traits, each having 10 chromosomes.

# GWAS SNPs
The number of RILs is 4,027 and the number of markers is 20,470,711.

Mixed linear model based association analysis.

### Output file format
`.mlma` (columns are chromosome, Markers, physical position, reference allele (the coded effect allele),
the other allele, frequency of the reference allele, SNP effect, standard error and p-value).

The significant QTLs are defined by -log10(p-value) > FDR (0.05).

```
library("IHW")
get_bh_threshold(GWAS.Results$p[!is.na(GWAS.Results$p)], 0.05)
```

For a few of traits, FDR threshold is not there, 5 is used.

360 output files are combinations of 36 traits and 10 chromosomes.
