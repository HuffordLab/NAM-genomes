#!/usr/bin/env perl
use strict;
use warnings;

#my $HQlib = $ARGV[0];
my $genome = "$ARGV[0].mod";
my $threads = 36;
my $path = $ARGV[1];

my $script_path = '~/las/git_bin/EDTA';
my $repeatmasker = '';
my $reclassify = "$script_path/util/classify_by_lib_RM.pl";
my $rename_by_list = "$script_path/util/rename_by_list.pl";
my $bed2gff = "$script_path/util/bed2gff.pl";

chdir $path;
        `perl $reclassify -seq $genome.EDTA.intact.fa -RM $genome.EDTA.intact.fa.out`;
        `perl $rename_by_list $genome.EDTA.intact.bed $genome.EDTA.intact.fa.rename.list 1 > $genome.EDTA.intact.bed.rename`;
        `mv $genome.EDTA.intact.bed.rename $genome.EDTA.intact.bed`;
        `perl $bed2gff $genome.EDTA.intact.bed`;
        `mv $genome.EDTA.intact.bed.gff $genome.EDTA.intact.gff`; #update intact.gff
        `cp $genome.EDTA.intact.gff ../`; #replace the intact gff that has no lib family info
