sfs="data/postcheck_out/post_sfs.txt"
param="data/postcheck_out/post_params.txt"

> $sfs
for i in `ls data/postcheck_out/`; do grep -h -B2 'SFS:' data/postcheck_out/$i | grep 'SFS:' | cut -d ':' -f2- | sed 's/^ //g' >> $sfs; done


cat <( ls data/postcheck_out/window* | sed 's/\_\_*/\t/g' | cut -f 2,4,6,8,10,12,14,16,18,20 | sort -u | sed 's/out\///g' ) <(ls data/postcheck_out/window* | sed 's/\_\_*/\t/g' | cut -f 3,5,7,9,11,13,15,17,19,21 ) > $param
