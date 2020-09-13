#!/usr/bin/env perl
use warnings;
use strict;
#get TE family percentage in genome
#Shujun Ou (shujun.ou.1@gmail.com)
#02/06/2020

my $usage = "perl get_TE_fam_pcnt.pl genome.RM.out.sum";
my $genome =$ARGV[0];
my $ID = $1 if $genome =~ /^([0-9a-z_]+)\..*/i;

open File, "grep % $genome|" or die $usage;
open Out1, ">$genome.fam" or die $usage;

my $go = 0;
while (<File>){
	($go = 1 and next) if /ID/;
	next unless $go == 1;
	s/^\s+//;
	my ($te, $cp, $bp, $pcnt) = (split);
	next unless $pcnt =~ /%/;
#	print "$te, $cp, $bp, $pcnt\n";
	print Out1 "$te\t$cp\t$bp\t$pcnt\t$ID\n";
	}
close File;
close Out1;

