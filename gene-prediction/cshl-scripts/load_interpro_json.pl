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
use JSON;

=head1 SYNOPSIS

load_interpro_json.pl  [options] interproscan.1.json interproscan.2.json
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species
    --debug
    --nowrite don't actually update the database

=head1 OPTIONS

=over 4

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name in the registry file

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=item B<--nowrite> 

   don't update the database

=back

=head1 ARGUMENTS

   interproscan json output files

=cut

my ($species, $registry);
my ($debug, $nowrite);
my $margin=undef;
{  							#Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
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

# get the interpro table
print STDERR "reading interpro table\n" if $debug;
my %iprids;
my $sql = "select interpro_ac,id from interpro";
my $sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($ipr,$id) = @$row;
  $iprids{$ipr}{$id} = 1;
}
$sth->finish;
# get analysis table
print STDERR "reading analysis table\n" if $debug;
my %analysisIdLUT;
$sql = "select analysis_id, logic_name from analysis";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($id,$logic_name) = @$row;
  print STDERR "analysis: $id, $logic_name\n" if $debug;
  $analysisIdLUT{uc($logic_name)} = $id;
}
$sth->finish;
# get translation lookup table
print STDERR "reading translation table\n" if $debug;
my %translationIdLUT;
$sql = "select translation_id, stable_id from translation";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
while (my $row = $sth->fetchrow_arrayref) {
  my ($tid,$stable_id) = @$row;
  $translationIdLUT{$stable_id} = $tid;
}
$sth->finish;

my $insert_interpro_sql = "insert into interpro (interpro_ac, id) VALUES (?,?)";
my @pf_cols = qw(translation_id seq_start seq_end hit_start hit_end hit_name analysis_id score evalue);
my $pf_cols_str = join(',',@pf_cols);
my $question_marks = join(',', map {'?'} @pf_cols);
my $insert_protein_feature_sql = "insert into protein_feature ($pf_cols_str), ($question_marks)";

my $insert_interpro_sth;
my $insert_protein_feature_sth;
 
unless ($nowrite){
  $insert_interpro_sth = $dbh->prepare($insert_interpro_sql) or die "cannot prepare $insert_interpro_sql\n";
  $insert_protein_feature_sth = $dbh->prepare($insert_protein_feature_sql) or die "cannot prepare $insert_protein_feature_sql\n";
}

for my $jsonfile (@ARGV) {
  print STDERR "Process $jsonfile\n" if $debug;
  open(my $fh, "<", $jsonfile);
  my $file_content = do {local $/; <$fh> };
  $file_content =~ s/}{/},{/g;
  my $ipr = decode_json $file_content;
  for my $res (@{$ipr->{results}}) {
    for my $xref (@{$res->{xref}}) {
      my $tid = $translationIdLUT{$xref->{id}};
      $tid or die "failed to get translation id for xref " . Dumper($xref);
      for my $match (@{$res->{matches}}) {
        my $hit_name = $match->{"model-ac"};
        my $ipr = $match->{signature}{entry}{accession};
        if (not exists $iprids{$ipr}{$hit_name}) {
          print STDERR "insert interpro $ipr $hit_name\n" if $debug;
          $insert_interpro_sth->execute($ipr, $hit_name) if $insert_interpro_sth;
        }
        my $score = $match->{score};
        my $analysis_id = $analysisIdLUT{$match->{signature}{signatureLibraryRelease}{library}};
        $analysis_id or die "failed to get analysis_id for signature " . Dumper($match->{signature});
        for my $location (@{$match->{locations}}) {
          my @pf = (
            $tid,
            $location->{start},
            $location->{end},
            $location->{hmmStart},
            $location->{hmmEnd},
            $hit_name,
            $analysis_id,
            $score,
            $location->{evalue}
          );
          print STDERR "insert protein_feature @pf\n" if $debug;
          $insert_protein_feature_sth->execute(@pf) if $insert_protein_feature_sth;
        }
      }
    }
  }
}


$insert_interpro_sth->finish if $insert_interpro_sth;
$insert_protein_feature_sth->finish if $insert_protein_feature_sth;
$dbh->disconnect;

  
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

