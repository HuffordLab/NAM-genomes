#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ($t2g_infile, $uc) = @ARGV;

my %t2g;
my %ss;
open (my $fh, "<", $t2g_infile);
while (<$fh>) {
  chomp;
  my ($t,$g) = split /\t/, $_;
  $t2g{$t}=$g;
  $ss{$g}{$t} = 1;
}
close $fh;

open ($fh, "<", $uc);
while (<$fh>) {
  chomp;
  my @x = split /\s+/, $_;
  next unless $x[0] eq 'H';
  my $tid = $x[8];
  $tid =~ s/_P/_T/;
  my $g = $t2g{$tid};
  $g or die "no t2g for tid='$tid'\n";
  delete $ss{$g} if exists $ss{$g};
}
close $fh;

for my $g (keys %ss) {
  for my $t (keys %{$ss{$g}}) {
    print "$t\t$g\n";
  }
}
