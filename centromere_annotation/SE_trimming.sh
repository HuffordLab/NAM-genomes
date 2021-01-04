#!/bin/bash

SE="$1"
wdir=$(dirname $SE)
cd $wdir
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
trim_galore --fastqc --gzip $SE -o $wdir
