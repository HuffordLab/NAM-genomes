#!/usr/local/bin/perl 

=head1 NAME

dump_transcripts.pl - Make fasta files of transcripts, and 3' and 5'
    regions, for all Genes or for Genes given as arguments.
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
#use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
#        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
#             ensembl-draw/modules ensembl-compara/modules );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;


#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Gene;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_transcripts.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --bylogicname	make a fasta file for each analysis logic_name
    --outfile 	the base filename
=cut

my ($species, $registry, $coding, $longest);
my (%exclude_gene,%exclude_analysispgm,%analysispgm,%exclude_clone,@logicNames,$outfile);
my $margin=undef;
{  #Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();
  my @exclude_analysispgm=();
  my @analysispgm=();
  my @exclude_clone=();
  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"bylogicname=s"=>\@logicNames
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"outfile=s"=> \$outfile
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

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 
#my $assembly_type=$ENS_DBA->get_MetaContainer()->get_default_assembly();

my @genes=@ARGV;
if (not @genes) {
  @genes=();
  for my $logicName (@logicNames) {
    my $gs = $gene_adaptor->fetch_all_by_logic_name($logicName);
    for my $g (@$gs) {
      push @genes, $g
    }
  }
}
# @genes or  @genes=@{$gene_adaptor->fetch_all_by_logic_name($bylogicname)};
print "retrieved gene count", scalar @genes;

#exit;

my %count;

my ($seqio_pr,%seqio,$seqio_cds,$seqio_cdna);

if ($outfile) {

	$seqio_pr= new Bio::SeqIO(-format => 'fasta',
                          -file => ">$outfile.evd.protein.fasta",
                         );
	 $seqio_cds= new Bio::SeqIO(-format => 'fasta',
                          -file => ">$outfile.evd.cds.fasta",
                         );
	 $seqio_cdna= new Bio::SeqIO(-format => 'fasta',
                          -file => ">$outfile.evd.cdna.fasta",
                         );
                          
}else{

	die "need outfile name";

}




foreach my $gene (@genes) {

    #print "! $geneid\n";

    $count{total_genes}++;
    if (@ARGV) {    # Stable ids GRMGnnnnnnn
        my $gid = $gene;
        eval { $gene = $gene_adaptor->fetch_by_stable_id($gid) };
        print STDERR "$@\n" and next if $@;
    }    #else {	#internal ids
         #  eval { $gene= $gene_adaptor->fetch_by_dbID($geneid); };
         #  print STDERR "gene_id $geneid:\n$@\n" and next if $@;
         #  # fails e.g. if gene has no exons
         #}

    next unless $gene;

    #my ($working_set_attrib) = @{$gene->get_all_Attributes('working-set')};
    #print "gene attrib is ", $working_set_attrib->code, "\n";
    #next if ($species =~ /zea|mays|maize/i && !$working_set_attrib);

    $count{qualified_genes}++;

    my @transcripts;

    eval { @transcripts = @{ $gene->get_all_Transcripts } };
    print STDERR "$@" && next if $@;

    foreach my $trans (@transcripts) {

        #print join "\t", ($trans->stable_id, $trans->spliced_seq, "\n");
        my $id = $trans->stable_id;

        my $cdna_seq = $trans->spliced_seq;

        my $seq_obj_cdna = Bio::Seq->new(
            -display_id => $id,
            -seq        => $cdna_seq,
        );

        $seqio_cdna->write_seq($seq_obj_cdna);
        if ( $trans->biotype eq 'protein_coding' ) {
            my $aa_obj = $trans->translate;

            eval { $seqio_pr->write_seq($aa_obj) };
            print STDERR "Cannot write seq for $id, $@" && next if $@;

            my $cdna_coding_start = $trans->cdna_coding_start;
            my $cdna_coding_end   = $trans->cdna_coding_end;
            my $seq_obj_cds       = Bio::Seq->new(
                -display_id => $id,
                -seq        => substr(
                    $cdna_seq,
                    $cdna_coding_start - 1,
                    $cdna_coding_end - $cdna_coding_start + 1
                )
            );

            $seqio_cds->write_seq($seq_obj_cds);
        }

        $count{qualified_transcripts}++;
    }


    
#  last;
}

for my $k (sort keys %count){
  print "$k = $count{$k}\n";
}



__END__


=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

