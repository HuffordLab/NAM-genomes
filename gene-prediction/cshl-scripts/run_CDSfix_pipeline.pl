#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;
use Pod::Usage;
use CDSFixer;

=head1 SYNOPSIS

run_CDSfix_pipeline.pl  [options] <list of gene stable identifiers>

 Options:
    --help     help message
    --man      full documentation
    --registry the registry file for database connections
    --species  which species to dump
    --prefix   prefix for gene annotation set, Zm00001eb
    --guide    guide file for conserved cds starts
    --debug    debug mode
    --nowrite  do not update database

=head1 OPTIONS

=over 4

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose transcripts are to be dumped

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=item B<--nowrite> 

   do not write to database

=item B<--prefix>

    prefix for gene ids e.g., Zm00001eb

=item B<--guide> 

   tab delimited file of upstream and internal conserved CDS start positions

=item B<--genes> 

   file listing one gene stable id per line

=back

=cut

my ($species, $registry, $debug, $prefix, $guide_file, $nowrite, $stage, $genes_file);

{ # Argument Processing
  my $help=0;
  my $man=0;

  GetOptions(
    "help|?"=>\$help
    ,"man"=>\$man
    ,"species=s"=>\$species
    ,"registry=s"=>\$registry
    ,"guide=s"=>\$guide_file
    ,"prefix=s"=>\$prefix
    ,"debug"=>\$debug
    ,"nowrite"=>\$nowrite
    ,"stage=i"=>\$stage
    ,"genes=s"=>\$genes_file
  ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
}

open (my $fh, "<", $genes_file);
my @genes;
for my $g (<$fh>) {
  chomp $g;
  push @genes, $g;
}
close $fh;

my %config = (
  registry  => $registry,
  species   => $species,
  genes     => \@genes,
  guideFile => $guide_file,
  prefix    => $prefix,
  debug     => $debug,
  nowrite   => $nowrite
);

my ($fixer,$err) = CDSFixer->new(%config);
$fixer or die $err;

# $fixer->runAll();

# or run specific steps
my @changes = ({},{},{},{},{},{});
$changes[1] = $fixer->findInternalStart({guided => 1}) if $stage == 1;
$changes[2] = $fixer->extendTranslation({guided => 1, prepend => 10, append => 10}) if $stage == 2;
$changes[3] = $fixer->replaceTranslation({minORF => 50}) if $stage == 3;
$changes[4] = $fixer->findInternalStart({MAX_EXON_INDEX => 2}) if $stage == 4;
$changes[5] = $fixer->extendTranslation({prepend => 10, append => 10, extend => 300}) if $stage == 5;

for my $i (1,2,3,4,5) {
  next unless $stage == $i;
  my %change = %{$changes[$i]};
  my $nGenes = keys %{$change{gene}};
  my $nTrans = keys %{$change{trans}};
  print "stage\t$i\tgenes\t$nGenes\ttrans\t$nTrans\n";
}
# tidy up db connections
$fixer->finish();

