## Steps to create PanGenome markers to use with the scripts

### Steps for processing PanGenome markers to use with ALLMAPS


File is located in CyVerse Data Commons that is accessible via iRods
```bash
icd /iplant/home/shared/panzea/genotypes/GBS/v27
iget Lu_2015_NatCommun_panGenomeAnchors20150219.txt.gz
```
or through direct link:

```
https://datacommons.cyverse.org/browse/iplant/home/shared/panzea/genotypes/GBS/v27/Lu_2015_NatCommun_panGenomeAnchors20150219.txt.gz
# wget/curl will not work
```

Once downloaded process the file downloaded as follows:

```bash
for f in {1..10}; do
 awk -v x=$f '$3==x && $7==0 {print $5"\tpg_"NR"\t"$8}' Lu_2015_NatCommun_panGenomeAnchors20150219.txt;
done > pb_anchors.txt
# bed file to extract sequence
for f in {1..10}; do
 awk -v x=$f '$3==x && $7==0 {print $5"\t"$6-50"\t"$6+50"\tpg_"NR"\t.\t+"}' Lu_2015_NatCommun_panGenomeAnchors20150219.txt;
done > pb_anchors.bed
# genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/Zea_mays.AGPv3.22.dna.genome.fa.gz
gunzip Zea_mays.AGPv3.22.dna.genome.fa.gz
# marker sequence
ml bedtools2
awk '$2>0' pb_anchors.bed > temp
mv temp pb_anchors.bed
bedtools getfasta -fi Zea_mays.AGPv3.22.dna.genome.fa -fo pb_anchors.fasta -bed pb_anchors.bed  -name
```

The final file `pb_anchors.txt` and `pb_anchors.fasta` will be used for AGP generation.


