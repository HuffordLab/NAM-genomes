#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $gff = shift @ARGV;

# parse the gff and output genes from the lookup tables
open (my $fh, "<", $gff);
while (<$fh>) {
  next if (/^#/);
  my @cols = split /\t/, $_;
  next unless $cols[2] eq 'mRNA';
  chomp $cols[8];
  my %attr;
  for my $kv (split /;/, $cols[8]) {
    my ($k, $v) = split /=/, $kv;
    $attr{$k} = $v;
  }
  print "$attr{ID}\t$attr{Parent}\n";
}
close $fh;

