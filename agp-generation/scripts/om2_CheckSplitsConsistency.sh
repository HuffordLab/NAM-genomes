#!/bin/bash
# checks the split onemap results for consistency
if [ "$#" -ne 1 ] ; then
echo "please provide split onemap result file";
exit 0;
fi
file="$1"
lg=$(cut -f 1 -d "," $file | sort | uniq)
markers=$(cut -f 2 -d "," $file | cut -f 1 -d "_" | sort | uniq -c | awk 'BEGIN{ORS=";"}{print $2"("$1")"}');
echo -e "${file%.*}\t${lg}\t${markers}";
# count number of markers in each linkage group and warn if they have less than 10
count=$(cat $file |wc -l)
if [ "${count}" -lt "10" ] ; then
echo "$(tput setaf 1)$file has less than 10 markers$(tput sgr0)";
fi
