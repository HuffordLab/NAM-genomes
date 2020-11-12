#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(devtools)
.libPaths("~/R/x86_64-pc-linux-gnu-library/3.5")
library(taxizedb)
library(phylostratr)
library(reshape2)
library(dplyr)
library(readr)
library(magrittr)
library(ggtree)
library(knitr)
weights=uniprot_weight_by_ref()
focal_taxid <- args[1]
strata <-
  uniprot_strata(focal_taxid, from=2) %>%
  strata_apply(f=diverse_subtree, n=5, weights=weights) %>%
  use_recommended_prokaryotes %>%
  add_taxa(c('4932', '9606')) %>%
  uniprot_fill_strata
strata_obj@data$faa[[focal_taxid]] <- '/ptmp/LAS/arnstrm/annotations-check/pyhylostrata/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.protein.fa'
