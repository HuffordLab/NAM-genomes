#!/usr/local/bin/perl 

=head1 NAME

set_canonical_transcripts_from_TRaCE.pl

TRaCE output format is a tab delimited file with columns
gene.stable_id
transcript.stable_id
rank
transcript length
CDS length
domain coverage (aa)
modified AED scores for each sample

This script only needs the first 3 columns. When rank is 1, set the canonical transcript id of the gene
To speed things up, fetch stable id to primary id lookup tables
=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

set_canonical_transcripts_from_TRaCE.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to dump
    --debug
    --nowrite don't actually update the database
    --trace output of TRaCE

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

=item B<--trace> 

   TRaCE output file

=back

=head1 ARGUMENTS


=cut

my ($species, $registry);
my ($debug, $nowrite, $traceFile);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
        ,"trace=s"=>\$traceFile
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  							#pod2usage(2) if $margin<0;
}
# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );
my $dbh = $ENS_DBA->dbc->db_handle;

# populate a lookup table of gene and transcript stable ids to gene/transcript id
my %IDLUT = (g => {}, t => {}, ct => {});
print STDERR "reading gene stable ids\n" if $debug;
my $sql = "select gene_id,stable_id,canonical_transcript_id from gene";
my $sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($gid,$stable_id,$ctid) = @$row;
  $IDLUT{g}{$stable_id} = $gid;
  $IDLUT{ct}{$stable_id} = $ctid;
}
$sth->finish;
print STDERR "reading transcript stable ids\n" if $debug;
$sql = "select transcript_id,stable_id from transcript";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($tid,$stable_id) = @$row;
  $IDLUT{t}{$stable_id} = $tid;
}
$sth->finish;

my $update_canonical_sql = "update gene set canonical_transcript_id=? where gene_id=?";

my $update_canonical_sth;
 
unless ($nowrite){
  $update_canonical_sth = $dbh->prepare($update_canonical_sql) or die "cannot prepare $update_canonical_sql\n";
}

print STDERR "reading TRaCE output\n" if $debug;
my $updates=0;
open (my $fh, "<", $traceFile);
while (<$fh>) {
  my ($g, $t, $rank, @etc) = split /\t/, $_;
  if ($rank == 1) {
    my $tid = $IDLUT{t}{$t};
    my $gid = $IDLUT{g}{$g};
    $tid or die "failed to get transcript id for $t";
    $gid or die "failed to get gene id for $g";
    my $ctid = $IDLUT{ct}{$g};
    if ($ctid != $tid) {
      print STDERR "setting canonical transcript for $g to $t\n" if $debug;
      $update_canonical_sth->execute($tid, $gid) if $update_canonical_sth;
      $updates++;
    }
  }
}
close $fh;

$update_canonical_sth->finish if $update_canonical_sth;
$dbh->disconnect;

print STDERR "set $updates canonical transcripts\n";
  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Andrew Olson <olson@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

