## Knob analyses

### Scripts

1. `Rep_content_dat_prep.R`: This script takes the Repeat_content_sum.csv output files for all lines
and combines them into a data set with more consistent columns. The output
dataset is called `NAM_array_coords.tsv`.

2. `Synteny_Knobs.R`: This script identifies knob synteny by comparing knob positions based on orthologous
genes. For the version for the main figure, only knobs that are syntenic to
B73 knobs that are >= 100kb are identified. For the supplemental version comparing
sequence-defined arrays and classical, cytological knobs, the orthologous gene sets
are limited to orthologs present in all lines. When any structural variants are identified
in the coordinates, the order in B73 is assumed to be the null, and the coordinates are
adjusted to match B73.

3. `NAM_plot_Supp.R_linux.sh` : This script modifies the gff files of repeats and TE's and subsetting to the elements that are within repeat arrays to reduce the required memory and time for further analysis. This is prep for NAM_plot_Supp.txt.

4. `NAM_plot_Supp.R` :This script generates the Supplement figure of the largest knobs with TE content. It takes the output from NAM_plot_Supp.R_linux and NAM_array_coords.tsv file as inputs.


### Data

* `NAM_array_coords.tsv` is the output from `Rep_content_dat_prep.R`.
* `NAM_array_coords_annotation_cyt_search_edit.csv` is the data from above with hand annotations checking for relationship with cytological knobs evident in FISH imagery.
* All scripts are in [`scripts`](scripts) folders and data in the [`assets`](assets) folder.
