<img src="/te-annotation/assets/Pan-EDTA_scheme.png" width=300 align=left>


#  De novo TE anno

## Identify raw TEs

```bash
type=tir #or ltr helitron
perl ~/las/git_bin/EDTA/EDTA_raw.pl --genome $genome --species Maize --type $type -t $threads
```
Generate de novo TE library for each genome and perform initial annotation

```bash
perl ~/las/git_bin/EDTA/EDTA.pl --genome $genome --species Maize -t $threads --cds $cds --curatedlib maizeTE02052020 --anno 1
```


## Make pan-genome TE library

Filter out single-copy annotations

```bash
# $genome.out is the RepeatMasker .out file generated from the last step, located in $genome.mod.EDTA.anno/
for i in `cat list.cds|awk '{print $1}'`; do
    perl ~/las/git_bin/EDTA/util/output_by_list.pl \
      1 \
     <(perl -nle 's/#.*//; print $_' $i.mod.EDTA.TElib.novel.fa) \
     1 \
     <(perl ~/las/git_bin/EDTA/util/find_flTE.pl $i.mod.out | \
     awk '{print $10}'| \
     sort| \
     uniq -c |\
     perl -nle 'my ($count, $id) = (split); if ($id=~/LTR/){next if $count<=2} else {next if $count ==1} print $_' |\
     awk '{print $2}') -FA > $i.mod.EDTA.TElib.novel.fa.real &
done
```

get classification info and convert #unknown to #DNA/Helitron

```bash
for j in *mod.EDTA.TElib.novel.fa; do 
    for i in `cat $j.real`; do 
        grep $i $j; 
    done| \
    perl -nle 's/#unknown/#DNA\/Helitron/; print $_' > $j.real.ori & 
done
```

aggregate novel TE libraries

```bash
i=0
for j in *real.ori; do
  i=$(($i+5000));
  perl ~/las/git_bin/EDTA/util/rename_TE.pl $j $i;
done > NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw
perl ~/las/git_bin/EDTA/util/rename_TE.pl NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw > NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw2
mv NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw2 NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw
```

remove redundant

```bash
nohup perl ~/las/git_bin/EDTA/util/cleanup_nested.pl \
    -in NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw \
    -cov 0.95 \
    -minlen 80 \
    -miniden 80 \
    -t 36 &
```
remove a number of false TEs and rename IDs

```
perl ~/las/git_bin/EDTA/util/output_by_list.pl \
     1 \
     NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw.cln \
     1 \
     rm.list2 \
     -ex \
     -FA > NAM.EDTA1.8.0.EDTA.TElib.novel.v2.fa.raw.cln2
perl ~/las/git_bin/EDTA/util/rename_TE.pl NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw.cln2 > NAM.EDTA1.8.0.EDTA.TElib.novel.fa
```
make comprehensive TE library

```
cat maizeTE02052020 NAM.EDTA1.8.0.EDTA.TElib.novel.v2.fa > NAM.EDTA1.8.0.MTEC02052020.TElib.fa
```


## Annotate pan-genome TE

Re-mask all genomes with the pan-genome TE lib

```bash
lib=NAM.EDTA1.8.0.MTEC02052020.TElib.fa
RepeatMasker -pa 36 -q -div 40 -lib $lib -cutoff 225 -gff $genome
for i in *fasta.out; do
 perl -nle 's/knob180_knob180/knob180/g;\
            s/CentC_CentC/CentC/g; \
            s/Unspecified/LTR\/CRM/g; \
            s/4\-12\-1_CL569186.1/CL569186.1/; \
            s/TR\-1_TR\-1/TR\-1/g;  \
            next if /TE_00018291_LTR/; \
            next if /TE_00013511/; \
            next if /TE_00016670_LTR/; \
            print $_' $i > $i.mod &
done
```

Re-run EDTA final on each genome

```bash
lib=NAM.EDTA1.8.0.MTEC02052020.TElib.fa
perl ~/las/git_bin/EDTA/EDTA.pl \
   --genome $genome \
   --species Maize \
   -t $threads \
   --step final \
   --anno 1 \
   --rmout $genome.out.mod \
   --curatedlib \$lib \
   --cds $cds
```


## calculate LAI

```bash
for i in `awk '{print $1}' list.cds`; do
 perl ~/las/git_bin/LTR_retriever/LAI \
    -genome $(echo $i|perl -nle 's/\..*//; print $_')/$i \
    -intact $(echo $i|perl -nle 's/\..*//; print $_')/$i.mod.EDTA.raw/LTR/$i.mod.pass.list \
    -all $i.out.mod \
    -q -iden 94.854 \
    -totLTR 76.34 \
    -t 2 &
done
```

## Summarize pan-genome TE

aggregate TE summary info

```bash
for i in *mod.EDTA.TEanno.sum; do
   cat <(echo $i| perl -nle 's/\..*//; print "$_\t${_}_cp\t${_}_bp\t${_}_pcnt"')\
       <(head -32 $i|grep -v -P "\-\-|=|total|masked"| perl -0777 -ne 's/\s+unknown/\nLTR_unknown/; print $_' |grep %) |\
   perl transpose3.pl -; \
done > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum
```

grep summary for each superfamily

```bash
cat head <(grep pcnt NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum) > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.pcnt.txt
cat head <(grep bp NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum) > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.bp.txt
cat head <(grep cp NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum) > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.cp.txt
```

extract family percent

```bash
for i in *mod.EDTA.TEanno.sum; do
  perl get_TE_fam_pcnt.pl $i &
done
```

aggregate into a big table

```bash
cat *mod.EDTA.TEanno.sum.fam | perl combine_TE_fam_pcnt.pl pcnt - > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam
cat *mod.EDTA.TEanno.sum.fam | perl combine_TE_fam_pcnt.pl bp - > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam.bp
perl -i -nle 's/#/_/g; print $_' NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam
perl -i -nle 's/#/_/g; print $_' NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam.bp
```

keep TE only

```bash
grep -v -P "CL569186.1|AF013103.1|\)n|cent|Cent|telo|knob|TR-1|osed|sela" NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam > NAM.EDTA1.9.0.MTEC02052020.TE.v1.1.anno.sum.fam
grep -v -P "CL569186.1|AF013103.1|\)n|cent|Cent|telo|knob|TR-1|osed|sela" NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.sum.fam.bp > NAM.EDTA1.9.0.MTEC02052020.TE.v1.1.anno.sum.fam.bp
```

count effective TEs

```bash
grep -v -P 'long_terminal_repeat|repeat_region|target_site_duplication' *mod.EDTA.TEanno.gff3 |\
perl -nle 'next unless s/ID=//; \
           my ($cla, $id)=(split)[2,8]; \
           $id=~s/.*; \
           Name=(.*);Classific.*/$1/; \
           $id=~s/;.*//; $id=~s/#/_/; \
           print "$id\t$cla"' |\
grep -v -P "\)n|A-rich|G-rich|begin|position" |\
sort -u > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.TEfam.list &
grep -v : NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.TEfam.list|wc -l
```

sum intact size for LTR TIR and Helitron in each genome

```bash
for j in `ls *fasta.mod.EDTA.intact.gff3`; do
   echo -n "$j ";
  for i in LTR_ /DT Heli; do
     grep ID $j | grep $i | awk '{print $1"\t"$4"\t"$5"\t"$3}'|\
     perl ~/las/git_bin/EDTA/util/combine_overlap.pl - |\
    perl ~/las/git_bin/EDTA/util/count_mask.pl -;
  done;
done |\
perl -ne 'chomp; print "\n" if /gff/; print "$_\t"' > NAM.EDTA1.9.0.MTEC02052020.intact.sum &
perl -i -nle 'next if /^$/; s/.*\///; s/\..*gff3//; print $_' NAM.EDTA1.9.0.MTEC02052020.intact.sum
cat <(echo "Genome LTR TIR Helitron") NAM.EDTA1.9.0.MTEC02052020.intact.sum > NAM.EDTA1.9.0.MTEC02052020.intact.sum.txt
````

get intact LTR info

```bash
for i in `ls *mod.EDTA.intact.gff3 |grep -v -P 'Ab10|AB10'`; do
  grep LTR_retrotransposon $i |\
  perl -nle 'my ($chr, undef, $supfam, $from, $to, undef, $str, undef, $info)=(split); \
             my $genome = $1 if $chr=~s/^(.*?)_//; \
             my ($id, $classification, $SO, $iden, $motif, $tsd)=($1, $2, $3, $4, $5, $6) if $info=~/Name=(.*);\
             Classification=(.*); \
             Sequence_ontology=(.*); \
             ltr_identity=(.*); \
             Method=structural; \
             motif=(.*);tsd=(.*)$/;\
            print "$genome\t$chr\t$supfam\t$classification\t$from\t$to\t$str\t$id\t$SO\t$motif\t$tsd\t$iden"';
done > NAM.26.intact.LTR.list &
```

plot figures in R

```
Rscript NAM_TE_summary.R
```


plot figures in R
```bash
Rscript Suppl. Fig. S4.R
Rscript Suppl. Fig. S5.R
Rscript Suppl. Fig. S6.R
```
