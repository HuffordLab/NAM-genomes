for f in /work/LAS/mhufford-lab/arnstrm/NAM_Canu1.8/chromosomes/*.fasta; do
targetdir="$(basename $f | cut -f 1 -d ".")/$(basename $f| cut -f 1 -d ".")_chr"
fasta=$(basename $f)
echo mkdir -p $targetdir
echo ln -s $f $targetdir/$fasta;
echo cd $targetdir
echo "echo /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/d-genetic-maps/break_and_index.sh $fasta > hisat.cmds";
echo "makeSLURMs.py 1 hisat.cmds";
echo cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/d-genetic-maps;
done
