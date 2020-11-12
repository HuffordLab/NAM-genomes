#!/usr/bin/perl
#
#Takes a ref file and a fasta file, and matches the masking in all the fasta entries according to the ref
#This has to be done chromosome-by-chromosome. sorry.
#
##Logan Kistler, Smithsonian Institution 2018. KistlerL@si.edu

use Getopt::Long;

GetOptions (
	'ref=s' => \$ref,
	'fasta=s' => \$fasta,
	'out=s' => \$out,
);

die "\nUsage: perl $0 --ref ref.fasta --fasta seqToMask.fasta --out out.fasta\n\n" unless ($ref and $fasta and $out);

open IN, $ref;
$line = <IN>;
$pos = 0;
while ($line = <IN>) {
	chomp $line;
	die "\n$line    looks like a fasta header, you're only allowed one of these\n\n" if ($line =~ /^>/);
	foreach $i (0..(length $line)-1) {
		$pos++;
		$nt = substr $line, $i, 1;
		$hold{$pos} = 1 if ($nt eq "N");
	}
}
close IN;

open IN, $fasta;
open OUT, ">$out";

while ($line = <IN>) {
	chomp $line ;
	if ($line =~ /^>/) {
		print OUT $line."\n";
		print "masking $line\n";
		$pos = 0;
	}
	else {
		foreach $i (0..(length $line)-1) {
			$pos++;
			if ($hold{$pos}==1) {
				substr $line, $i, 1, "N";
			}
		}
		print OUT $line."\n";
	}
}
