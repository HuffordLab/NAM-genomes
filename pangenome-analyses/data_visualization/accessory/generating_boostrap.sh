for k in {1..10}; do for j in {1..100}; do for i in `ls *list.txt|shuf`; do cat $i >> temp.$j.$k; sort -u temp.$j.$k|wc -l; 
done |bash transpose.sh - > result.$j.$k; rm temp.$j.$k; done & done
