#! /usr/bin/perl -w

# pull_out_caononical_transcript_coded_from_gff.pl
# 28 April 2020
# Candy Hirsch

use strict;
use Getopt::Std;

my $usage = "\n$0 -i input_gff_file -l input_list_of_canonical_transcripts -o output_gff_file\n";

our ($opt_i, $opt_l, $opt_o, $opt_h);
getopts("i:l:k:o:h") or die "$usage";

if ( (!(defined $opt_i)) || (!(defined $opt_l)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
  exit;
}

open (my $GFFIn_fh, '<', $opt_i) || die "\nCannot open $opt_i\n\n";
open (my $list_fh, '<', $opt_l) || die "\nCannot open $opt_l\n\n";
open (my $GFFOut_fh, '>', $opt_o) || die "\nCannot open $opt_o\n\n";

my %canonical;
while (my $line = <$list_fh>) {
  chomp $line;
  $canonical{$line} = 0;
}

while (my $line = <$GFFIn_fh>) {
  chomp $line;

  if ($line =~ /^#/) {
    next;
  }

  my (undef, undef, $type, undef, undef, undef, undef, undef, $info) = split ("\t", $line);
  my ($name, undef) = split (";", $info);
  $name =~ s/ID=//;

  if ($type eq 'mRNA') {
    
    if (defined $canonical{$name}) {
      ++$canonical{$name};
      print $GFFOut_fh "$line\n";
    }
  }
}


foreach my $keys (keys %canonical) {
  if ($canonical{$keys} != 1) {
    print "$keys\t$canonical{$keys}\n";
  }
}


close $GFFIn_fh;
close $GFFOut_fh;
close $list_fh;
