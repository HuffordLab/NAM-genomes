# NLR-genes

## Methods

The NLRs were extracted from the nucleotide genomic sequences using `NLR-Annotator` (Steuernagel et al. 2018) and from proteomes using `hmmalign` with reference HMM of the grass NB-ARC (Bailey et al. 2018).
Additionally, NLRs and NLR-IDs were characterized in the _Brachypodium_ (Gordon et al. 2017) and Maize annotations using the [plant_rgenes pipeline](https://github.com/krasileva-group/plant_rgenes) (Sarris et al. 2016) (`e-value` cutoff `1e-3`).
The number of NB-ARC containing proteins was compared to those previously identified in Arabidopsis (Van et al. 2019) and plotted using `R` package `ggplot2` (Wickham 2016).
The NB-ARC domain  alignment was manually curated for the presence of NB-ARC domain functional motifs including `Walker A`, `WALKER-B`, `RNBS-C`, `GLPL` and `RNBS-D`.
The NLR phylogeny was determined using `RAxML` MPI (v8.2.9, `-f a`, `-x 12345`, `-p 12345`, `-# 100`, `-m PROTCATJTT`).
The phylogeny was visualized and re-rooted on the longest internal branch in iTOL (Letunic and Bork 2016).    

## Other Data:

For the ITOL scripts you need to append to start the output the formating that corresponds to the iTOL annotation pages which can be found here: https://itol.embl.de/help.cgi
Link to trees shared: https://itol.embl.de/shared/xCJbI9ndshEK
