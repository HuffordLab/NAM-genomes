ABC based analysis of Maize demography and selection. Modifying priors, number of simulations, etc. is all done by modifying the Snakefile. Run with `snakemake -j X`, where X is the desired number of parallel jobs to run. For HPC, `bash submit.sh` can be used.

Once all simulations are completed, the final files for analysis can be generated with
`bash src/make_df.sh`, which will produce  `./pop.txt` and `./pop_sfs.txt`, which are used as input for ABC. 


