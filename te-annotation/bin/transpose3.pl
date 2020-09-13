#!usr/bin/perl -w
use strict;
#transpose a table. Memory intensive! 2GB size table consumes 50+GB memory for 10 min.
#Usage: perl transpose3.pl table >tr_table
#Author: Shujun Ou 2014-02-14
#Version: v1.0

my $i=0;
my @ori;
my @trans;
while (<>){
	chomp;
	@{$ori[$i]}=(split);
	for (my $j=0; $j<=$#{$ori[$i]}; $j++){
		$trans[$j][$i]=$ori[$i][$j];
		}
	@{$ori[$i]}='';
	$i++;
	}

for (my $a=0; $a<=$#trans; $a++){
	for (my $b=0; $b<=$#{$trans[$a]}; $b++){
		print "$trans[$a][$b]\t";
		}
	print "\n";
	@{$trans[$a]}='';
	}


