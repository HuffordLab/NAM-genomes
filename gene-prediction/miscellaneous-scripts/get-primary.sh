#!/bin/bash
gff=$1
grep "canonical_transcript=1" $gff |cut -f 9 | cut -f 1-2 -d ";" | sed 's/;Parent=gene:/\t/g' |sed 's/^ID=transcript://g' |awk '{print $2"\t"$1}'
