#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ($prefix, $uc, $t2g, $protein_fasta, $cds_fasta, $gff) = @ARGV;

my %transcripts;
# gather transcripts that have a hit from uclust output file
open (my $fh, "<", $uc);
while (<$fh>) {
  chomp;
  my @cols = split /\t/, $_;
  if ($cols[0] eq 'H') {
    $cols[8] =~ s/\sgene.*//;
    $cols[8] =~ s/_P/_T/;
    $cols[8] =~ s/\s*$//;
    $transcripts{$cols[8]} = 1;
  }
}
close $fh;

# read the transcript to gene lookup table (already made from gff)
my %keep = (t => {}, g=> {});
my %g2t;
open ($fh, "<", $t2g);
while (<$fh>) {
  chomp;
  my ($t, $g) = split /\t/, $_;
  $g2t{$g}{$t}=1;
  $keep{g}{$g} = 1 if $transcripts{$t};
}
close $fh;
for my $g (keys %{$keep{g}}) {
  for my $t (keys %{$g2t{$g}}) {
    $keep{t}{$t}=1;
  }
}

# parse the gff and output genes from the lookup tables
open ($fh, "<", $gff);
open (my $outfh, ">", "$prefix.gff");
my $prev = 'not likely to be the first line in the gff';
while (<$fh>) {
  if (/^#/) {
    print $outfh $_ if $_ ne $prev;
    $prev = $_;
    next;
  }
  my @cols = split /\t/, $_;
  chomp $cols[8];
  my %attr;
  for my $kv (split /;/, $cols[8]) {
    my ($k, $v) = split /=/, $kv;
    $attr{$k} = $v;
    if ($k eq 'Parent') {
      $attr{$k} = {};
      for my $p (split /,/, $v) {
        $attr{$k}{$p}=1;
      }
    }
  }
  if (($cols[2] eq 'gene' and $keep{g}{$attr{ID}})
   or ($cols[2] eq 'mRNA' and $keep{t}{$attr{ID}})
   or (exists $attr{Parent} and hasOne($keep{t},keys %{$attr{Parent}}))
  ) {
    print $outfh $_;
    $prev = $_;
  }
}
close $outfh;
close $fh;

# filter the fasta files
open ($fh, "<", $protein_fasta);
open ($outfh, ">", "$prefix.protein.fasta");
my $id;
while (<$fh>) {
  if (/^>(\S+)/) {
    $id = $1;
    $id =~ s/_P/_T/;
  }
  print $outfh $_ if ($keep{t}{$id});
}
close $outfh;
close $fh;

open ($fh, "<", $cds_fasta);
open ($outfh, ">", "$prefix.cds.fasta");
while (<$fh>) {
  if (/^>(\S+)/) {
    $id = $1;
  }
  print $outfh $_ if ($keep{t}{$id});
}
close $outfh;
close $fh;

sub hasOne {
  my ($hsh, @keys) = @_;
  for my $k (@keys) {
    return 1 if exists $hsh->{$k};
  }
  return 0;
}