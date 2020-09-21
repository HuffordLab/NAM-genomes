maf_array=($( ls -d ./analyses/last/net_axt/net_maf/*1to1.maf ))
combined_maf=./analyses/last/net_axt/net_maf/combined.maf

cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp

multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf

for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done

# and filter mafs so all blocks have Zea mays and are at least 20 bp long
mafFilter -minCol=20 -needComp="$ref_name" $combined_maf > "$combined_maf".filtered
