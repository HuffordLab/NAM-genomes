#!/bin/bash
nam=$1
rm -rf tsp_work combined.* weights.txt unmapped chr*pdf
sort ${nam}_GG-mapped.csv |uniq > temp
mv temp ${nam}_GG-mapped.csv
../../gg4_clean-markers.sh ${nam}_GG-mapped.csv
