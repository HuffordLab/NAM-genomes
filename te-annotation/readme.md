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
cat maizeTE02052020 NAM.EDTA1.8.0.EDTA.TElib.novel.v2.fa > NAM.EDTA1.9.0.MTEC02052020.TElib.fa
```


## Annotate pan-genome TE

Re-mask all genomes with the pan-genome TE lib

```bash
lib=NAM.EDTA1.9.0.MTEC02052020.TElib.fa
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
lib=NAM.EDTA1.9.0.MTEC02052020.TElib.fa
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

plot figures in R
```bash
Rscript Suppl. Fig. S4.R
Rscript Suppl. Fig. S5.R
Rscript Suppl. Fig. S6.R
```
