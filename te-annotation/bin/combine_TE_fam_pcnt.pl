#!/usr/bin/env perl
use warnings;
use strict;
#Usage: compile family info of different lines into a big table
#Shujun Ou (shujun.ou.1@gmail.com)
#02/06/2020

my $usage = "cat *.sum.fam | perl combine_TE_fam_pcnt.pl [bp|pcnt] -";
my $type = $ARGV[0];
$type = 'bp' unless defined $type;
my %summary;
open IN, "<$ARGV[1]" or die $usage;

while (<IN>){
	my ($te, $cp, $bp, $pcnt, $id) = (split);
	next unless defined $id;
	if ($type eq 'pcnt'){
		${$summary{$te}}{$id} = $pcnt;
		}
	elsif ($type eq 'bp'){
		${$summary{$te}}{$id} = $bp;
		}
	}
close IN;

	
print "TE_fam";
#use "CentC" to extract genome id info
foreach my $id (sort {$a cmp $b} (keys %{$summary{"CentC"}})){
	print "\t$id";
	}
print "\n";

foreach my $te (keys %summary){
	print "$te";
	foreach my $id (sort {$a cmp $b} (keys %{$summary{"CentC"}})){
		my $pcnt = 0;
		$pcnt = ${$summary{$te}}{$id} if exists ${$summary{$te}}{$id};
		$pcnt = $pcnt/100 if $pcnt =~ s/%//; #convert percent to decimal
		print "\t$pcnt";
		}
	print "\n";
	}
