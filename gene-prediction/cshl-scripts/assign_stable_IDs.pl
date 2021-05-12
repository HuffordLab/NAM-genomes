#!/bin/env perl

=head1 NAME

assign_stable_IDs.pl

Assigns stable gene identifiers of the form <prefix><5 digits>0
The 5 digits are a number based on the order of the genes along the chromosomes

The prefix is supplied as an argument. So, for example: Zm00001eb000030 means:
  the third gene (in order) 
  of the 2nd annotation release (b)
  on the 5th version of the assembly (e)
  of the B73 genotype (00001)
  of the species Zea mays (Zm)

=cut

use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use POSIX;
use Scalar::Util qw (reftype);

#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

use Getopt::Long;
use Pod::Usage;

=head1 SYNOPSIS

assign_stable_IDs.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to dump
    --prefix      prefix for gene annotation set, Zm00001eb
    --bylogicname	only work on genes udner these analysis logic_name(s)
    --debug             debug mode
    --nowrite     do not update database

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

=item B<--bylogicname> 

   comma separated list of analysis logic_names to assign IDs to

=back

=cut

my ($species, $registry, $debug, $prefix, @logicNames, $nowrite);

{  							#Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
        ,"bylogicname=s"=>\@logicNames
        ,"prefix=s"=>\$prefix
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

my @tables = qw(gene exon transcript translation);
my %has_description = (
  gene => 1,
  transcript => 1
);
# prepare update statements
my %update_sth;
for my $table (@tables) {
  my $extra = $has_description{$table} ? ', description=?' : '';
  my $sql = "update $table set stable_id=? $extra where ${table}_id=?";
  $update_sth{$table} = $dbh->prepare($sql) or die "cannot prepare SQL '$sql'\n" . $dbh->errstr;
}

# get genes in genome order
my $logicNameStr = join(",", map { "'$_'" } @logicNames);
my $sql = qq{
select g.gene_id,g.stable_id, sr.name, g.seq_region_start, g.seq_region_end
from gene g, analysis a, seq_region sr
where g.analysis_id = a.analysis_id and a.logic_name in ($logicNameStr)
and g.seq_region_id = sr.seq_region_id
order by g.seq_region_id, g.seq_region_start};

print STDERR "gene sql : $sql\n" if $debug;

my $sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
my @genes;
my $gene_idx=0;
while (my $row = $sth->fetchrow_arrayref) {
  $gene_idx++;
  my ($gid, $old_sid, $sr, $st, $e) = @$row;
  my $new_sid = sprintf("%s%05d0",$prefix,$gene_idx);
  push @genes, [$gid, $old_sid, $new_sid, $sr, $st, $e];
}
$sth->finish;

# fetch transcripts
$sql = "select transcript_id, gene_id, seq_region_start, seq_region_end, seq_region_strand from transcript";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
my %g2t;
while (my $row = $sth->fetchrow_arrayref) {
  my ($tid, $gid, $st, $e, $ori) = @$row;
  $g2t{$gid}{$tid} = [$st, $e, $ori];
}
$sth->finish;

# fetch translations
$sql = "select translation_id, transcript_id from translation";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
my %tr2tl;
while (my $row = $sth->fetchrow_arrayref) {
  my ($tl, $tr) = @$row;
  $tr2tl{$tr} = $tl;
}
$sth->finish;

# fetch exon_transcript table
$sql = "select exon_id, transcript_id, rank from exon_transcript";
$sth = $dbh->prepare($sql) or die "Error:" . $dbh->errstr . "\n";
$sth->execute or die "Error:" . $sth->errstr . "\n";
my %t2e;
while (my $row = $sth->fetchrow_arrayref) {
  my ($eid, $tid, $rank) = @$row;
  $t2e{$tid}{$eid} = $rank;
}
$sth->finish;

# iterate over the genes and update stable ids of gene and associated transcripts, translations, and exons
for my $gene (@genes) {
  print join("\t", @$gene),"\n";
  my ($gid, $old_stable_id, $gene_stable_id, @etc) = @$gene;
  $update_sth{gene}->execute($gene_stable_id, $gene_stable_id, $gid) unless $nowrite;
  print STDERR "update gene set stable_id='$gene_stable_id' where gene_id=$gid\n" if $debug;
  my $trHsh = $g2t{$gid};
  my @tids = sort {
    ($trHsh->{$a}[2] == -1
    and $trHsh->{$b}[1] <=> $trHsh->{$a}[1] or $trHsh->{$a}[0] <=> $trHsh->{$b}[0])
    or ($trHsh->{$a}[0] <=> $trHsh->{$b}[0] or $trHsh->{$b}[1] <=> $trHsh->{$a}[1])
  } keys %$trHsh;
  my $tidx = 0;
  my %seen; # exons we have already picked a stable id for
  for my $tid (@tids) {
    $tidx++;
    my $tr_stable_id = sprintf("%s_T%03d", $gene_stable_id, $tidx);
    $update_sth{transcript}->execute($tr_stable_id, $tr_stable_id, $tid) unless $nowrite;
    print STDERR "update transcript set stable_id='$tr_stable_id' where transcript_id=$tid\n" if $debug;
    if ($tr2tl{$tid}) {
      my $tl_stable_id = $tr_stable_id;
      $tl_stable_id =~ s/_T/_P/;
      $update_sth{translation}->execute($tl_stable_id, $tr2tl{$tid}) unless $nowrite;
      print STDERR "update translation set stable_id='$tl_stable_id' where translation_id=$tr2tl{$tid}\n" if $debug;
    }
    my @new_exons_by_rank = sort {$t2e{$tid}{$a} <=> $t2e{$tid}{$b}} grep {not $seen{$_}} keys %{$t2e{$tid}};
    for my $eid (@new_exons_by_rank) {
      my $exon_stable_id = "$tr_stable_id.exon$t2e{$tid}{$eid}";
      $seen{$eid} = $exon_stable_id;
      $update_sth{exon}->execute($exon_stable_id, $eid) unless $nowrite;
      print STDERR "update exon set stable_id='$exon_stable_id' where exon_id=$eid\n" if $debug;
    }
  }
}

for my $table (@tables) {
  $update_sth{$table}->finish;
}
