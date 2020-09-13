#!/usr/bin/env perl
use strict;
use warnings;
# this script is to generate slurm scripts for batch analyses
# uncomment the necessary lines for different purpose
# Usage: perl this_script.pl genome_list
# Author: Shujun Ou (Sept. 2020) (shujun.ou.1@gmail.com)

while (<>) {
	chomp;
	next if /^$/;
	my ($info, $cds) = (split);
	my ($path, $genome)=($1, $2) if $info=~/(.*)\/(.*)/;
	$genome=$info;
	my $ID = $1 if $info =~ /^([0-9a-z_]+)\./i;

	my $type = "ltr";
#	my $type = "tir";
#	my $type = "helitron";
#	my $type = "EDTA";
#	my $type = "RM";
#	my $type = "final";

	my $jobtype="$ID.$type";
	my $threads=36;
	open Qsub, ">worker.$jobtype.qsub";
	print Qsub "#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 4:00:00
##SBATCH -t 1-0:00:00
#SBATCH --ntasks-per-node $threads
#SBATCH --mem=300GB

conda activate EDTA

workdir=/home/MaizeNAM/TE_anno

mkdir $ID
cd \$workdir/$ID

ln -s ~/jfw/TE/MaizeNAM/NAM_canu1.8/verified/$genome $genome
ln -s /home/oushujun/jfw/TE/MaizeNAM/TE_anno/$cds $cds

# de novo library
#perl ~/las/git_bin/EDTA/EDTA_raw.pl --genome $genome --species Maize --type $type -t $threads
#perl ~/las/git_bin/EDTA/EDTA.pl --genome $genome --species Maize -t $threads --cds $cds --curatedlib maizeTE02052020

# re-mask genomes with the pan-genome lib
lib=NAM.EDTA1.9.0.MTEC02052020.TElib.fa
#ln -s /work/LAS/jfw-lab/oushujun/TE/MaizeNAM/TE_anno/\$lib \$lib
#RepeatMasker -pa 36 -q -div 40 -lib \$lib -cutoff 225 -gff $genome

# re-run EDTA final
lib=NAM.EDTA1.9.0.MTEC02052020.TElib.fa
#/usr/bin/time -v perl ~/las/git_bin/EDTA/EDTA.pl --genome $genome --species Maize -t $threads --step final --anno 1 --rmout /home/oushujun/jfw/TE/MaizeNAM/TE_anno/NAM.EDTA1.8.0/pan_TE_lib_mask/$genome.out --curatedlib \$lib --cds $cds


";
close Qsub;
	}
