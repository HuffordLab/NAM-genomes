package CDSFixer;

=pod

=head1 CDSFixer

Library to fix CDS issues in ensembl cores

=head1 Example usage

my %config = (
  registry => '/path/to/reg.pm',
  species  => 'species identifer in reg.pm',
  genes    => \@list_of_gene_ids,
  guide    => '/path/to/guided_fixes.txt',     # optional tab delimited file (upstream|internal), translationStableID, new CDS start position
  debug    => 1                                # set debug level
);

my ($fixer,$err) = CDSFixer->new(%config);
$fixer or die "$err";

$fixer->runAll();

# or run specific steps
$fixer->findInternalStart({guided => 1});
$fixer->extendTranslation({guided => 1});
$fixer->replaceTranslation();
$fixer->findInternalStart();
$fixer->extendTranslation();

# tidy up db connections
$fixer->finish();

=cut

use strict;
use warnings;
use Data::Dumper;

use POSIX;

use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

=pod

=head1 METHODS

=item new

=over

Constructor - initialize db connection, fetch some data needed in other steps

=cut

sub new {
  my $class = shift;
  my %init = @_;

  my $self = bless {}, $class;

  return $self->initialize(%init);
}

=pod

=item initialize

Connects to db, gets genes and fetches exon transcript info for requested genes.
Called automatically by the constructor, must be explicitly called if you run more than once from the same CDSFixer object

=cut

sub initialize {
  my $self = shift;
	my %init = @_;
  for my $param ('registry','species','genes') {
    $init{$param} or return (undef, error("'$param' parameter is not defined"));
  }
  print STDERR "debug = '$init{debug}'\n";
  $init{debug} ||= 0;
  for my $param (keys %init) {
    $self->{$param} = $init{$param};
  }
  if ($init{guideFile}) {
    print STDERR "will read guid from $init{guideFile}\n";
    return (undef, error("'guide' file does not exists '$init{guideFile}")) unless -e $init{guideFile};
    open (my $guide_fh, "<", $init{guideFile});
    while (<$guide_fh>) {
      chomp;
      my ($direction, $pid, $position, @etc) = split /\t/, $_;
      $self->{guide}{$direction}{$pid} = $position;
    }
    close $guide_fh;
  }

  # Load the ensembl registry file and connect to database
  Bio::EnsEMBL::Registry->load_all( $init{registry} );
  my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $init{species}, 'core' );
  $ENS_DBA or return (undef, error("No core DB for species '$init{species}' in registry $init{registry}"));

  my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
	
  @{$self->{genes}} = map { $gene_adaptor->fetch_by_stable_id($_) } @{$init{genes}};

  print STDERR "#INFO will process ",scalar @{$self->{genes}}," genes\n" if $self->{debug};

  my $dbh = $ENS_DBA->dbc->db_handle;

  # get the exon transcript phase for the genes we will be working on
  my %exon_transcript_phase; # {exon}{transcript} = phase
  my $geneIdList = join(',', map {$_->dbID} @{$self->{genes}});
  my $sth = $dbh->prepare(qq{
    select et.exon_id, et.transcript_id, et.rank, e.phase from exon_transcript et, exon e, transcript t
    where et.exon_id = e.exon_id and et.transcript_id = t.transcript_id and t.gene_id IN ($geneIdList)
    order by et.transcript_id, et.rank DESC
  });
  $sth->execute();
  $self->{exonCheck} = {};
  my %lastExonLUT;
  
  while (my $row = $sth->fetchrow_arrayref) {
    my ($e,$t,$rank,$phase) = @$row;
    if (exists $lastExonLUT{$t}) {
      $self->{exonCheck}{isNotLast}{$e} = $t;
    }
    else {
      $lastExonLUT{$t}=$e;
    }
    if ($rank > 1) {
      $self->{exonCheck}{isNotFirst}{$e} = $t;
    }
    $exon_transcript_phase{$e}{$t} = $phase;
  }
  $sth->finish;
  print STDERR "read exon transcript phase\n" if $self->{debug};
  $self->{exon_transcript_phase} = \%exon_transcript_phase;

  # prepare SQL statements
  $self->{sth} = prepareSQL($dbh);
  $self->{dbh} = $dbh;
  return $self;
}

sub finish {
  my $self = shift;
  for my $sth (values %{$self->{sth}}) {
    $sth->finish;
  }
  $self->{dbh}->disconnect;
}

sub error {
  my $message = shift;
  return "CDSFixer: ERROR: $message\n";
}


sub prepareSQL {
  my $dbh = shift;
  my @exon_fields = qw(exon_id seq_region_id seq_region_start seq_region_end seq_region_strand phase end_phase is_current is_constitutive stable_id version);
  my $exon_fields_string = join(',',@exon_fields);
  my $question_marks = join(',', map {'?'} @exon_fields);

  my %sql = (
    select_exon               => "select $exon_fields_string from exon where exon_id=?",
    insert_exon               => "insert into exon ($exon_fields_string) VALUES ($question_marks)",
    update_exon_phase         => "update exon set phase=? where exon_id=?",
    update_exon_transcript    => "update exon_transcript set exon_id=? where exon_id=? and transcript_id=?",
    update_translation        => "update translation set start_exon_id=?, seq_start=?, end_exon_id=? where translation_id=?",
    update_translation_start  => "update translation set start_exon_id=?, seq_start=? where translation_id=?",
    update_translation_end    => "update translation set end_exon_id=?, seq_end=? where translation_id=?",
    update_translation_seq    => "update translation set seq_start=?, seq_end=? where translation_id=?",
    update_gene_bounds        => "update gene set seq_region_start=?, seq_region_end=? where gene_id=?",
    update_transcript_bounds  => "update transcript set seq_region_start=?, seq_region_end=? where transcript_id=?",
    update_exon_bounds        => "update exon set seq_region_start=?, seq_region_end=? where exon_id=?",
    delete_translation        => "delete from translation where transcript_id=?",
    make_transcript_noncoding => "update transcript set biotype='non_coding' where transcript_id=?",
    make_gene_noncoding       => "update gene set biotype='non_coding' where gene_id=?"
  );

  my %sth; 
  for my $statement (keys %sql) {
    $sth{$statement} = $dbh->prepare($sql{$statement}) or die error("cannot prepare $sql{$statement}");
  }
  
  return \%sth;
}

sub findInternalStart {
  my ($self,$params) = @_;
  my $MIN_TRANSLATION_LENGTH = $params->{MIN_TRANSLATION_LENGTH} || 50;
  my $MAX_EXON_INDEX = $params->{MAX_EXON_INDEX} || 99;
  my %changes;
  foreach my $gene (@{$self->{genes}}) {
    my %updates;
    next unless $gene->biotype eq 'protein_coding';
    my $strand = $gene->strand;
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      next unless $transcript->biotype eq 'protein_coding';
      my $translation = $transcript->translation;
      next unless $translation;
      my $translation_id = $translation->dbID;
      my $idxm=0;
      my $aa;
      if ($params->{guided}) {
        next unless $self->{guide}{internal}{$translation->stable_id};
        $idxm = $self->{guide}{internal}{$translation->stable_id};
        $aa = $transcript->translate->seq;
      }
      else {
        $aa = $transcript->translate->seq;
        $idxm = index( uc($aa), 'M', 0);
        next unless ($idxm > 0); # skip translations that start with M or don't contian an M 
      }
      $idxm += 1;
      next if (length($aa) - $idxm < $MIN_TRANSLATION_LENGTH);

      my $trmapper = $transcript->get_TranscriptMapper;
      my $endExonID = $translation->end_Exon->dbID;

      my @genomic_coords = $trmapper->pep2genomic( $idxm, $idxm );
      my $Met_start_genomic;
      my $start_exon_start_phase;

      if( scalar @genomic_coords == 0 ){
        warn("No genomic coord found for M at $idxm, skip\n");
        next;
      }

      if ($strand > 0){
        my @starts = sort map {$_->start} @genomic_coords;
        $Met_start_genomic = shift @starts;
      } else {
        my @ends = reverse sort map {$_->end} @genomic_coords;
        $Met_start_genomic = shift @ends;
      }

      my @fiveUTRexonIDs2update;
      my ($met_start_ExonID, $start_exon_start, $exon_start_phase);
      my @ordered_Exons = $strand > 0
      ? @{$transcript->get_all_Exons}
      : sort {$b->seq_region_start <=> $a->seq_region_start} @{$transcript->get_all_Exons};
      my $Exon;
      my $exon_idx=0;
      while ( $Exon = shift @ordered_Exons){
        $exon_idx++;
        my $exon_gstart = $Exon->seq_region_start;
        my $exon_gend = $Exon->seq_region_end;
        $exon_start_phase = $Exon->phase;
        my $exon_seq = $Exon->seq->seq;

        if ($Met_start_genomic >= $exon_gstart && $Met_start_genomic <= $exon_gend) {
          $met_start_ExonID = $Exon->dbID;
          $start_exon_start = $strand > 0 ? $Met_start_genomic-$exon_gstart+1:$exon_gend-$Met_start_genomic+1;
          if ($strand > 0){
            $start_exon_start_phase = $Met_start_genomic == $exon_gstart ? 0 : -1;
          }else{
            $start_exon_start_phase = $Met_start_genomic == $exon_gend ? 0 : -1;
          }
          last;
        }

        push @fiveUTRexonIDs2update, $Exon->dbID if $strand>0;
      }
      @fiveUTRexonIDs2update = map{ $_->dbID } @ordered_Exons if $strand < 0;

      unless( $met_start_ExonID && defined $start_exon_start_phase){
        warn("ERROR: no valid exon found and start found for genomic coord   $Met_start_genomic\n");
        next;
      }

      # print STDERR "$translation_stable_id ($translation_id), old start $translation_old_start, startExonID $met_start_ExonID, Met start in startExon $start_exon_start, startPhase ($exon_start_phase -> $start_exon_start_phase) exon_idx $exon_idx\n";
      next if ($exon_idx > $MAX_EXON_INDEX);
      $changes{gene}{$gene->dbID}=1;
      $changes{trans}{$transcript->dbID}=1;
      unless($self->{nowrite}) {
        $updates{$met_start_ExonID}{$transcript->dbID} = {
          startExon => $met_start_ExonID, # this could change if a new exon needs to be created
          phase => $start_exon_start_phase,
          seqStart => $start_exon_start,
          endExon => $endExonID, # this could change if a new exon needs to be created
          translationID => $translation_id
        };
      }
    }
    $self->doUpdates(\%updates);
  }
  return \%changes;
}

sub extendTranslation {
  my ($self,$params) = @_;
  my $extendM = $params->{extendM};
  my $append = $params->{append} || 0;
  my $prepend = $params->{prepend} || 0;
  my $extend = $params->{extend} || 2001;
  my %guide = %{$self->{guide}{upstream}};
  my $sth = $self->{sth};
  my %changes;
  for my $gene (@{$self->{genes}}) {
    next unless $gene->biotype eq 'protein_coding';
    my $geneID = $gene->dbID;
    my $geneStart = $gene->start;
    my $geneEnd = $gene->end;
    my $updateGene=0;
    my @transcripts = @{$gene->get_all_Transcripts()};
    my (%adjustedStart, %adjustedEnd);
    for my $trans (@transcripts) {
      my $transID = $trans->dbID;
      next unless $trans->biotype eq 'protein_coding';
      my $translation = $trans->translation;
      next unless $translation;
      my $translationID = $translation->dbID;
      next if ($params->{guided} and not $guide{$translation->stable_id});
      my $firstAA = substr($translation->seq,0,1);
      my @exons = @{$trans->get_all_Exons()};
      my $startExonID = $trans->start_Exon->dbID;
      my $updateFirstExon=0;
      my $updateLastExon=0;
      if (not exists $adjustedStart{$startExonID}) {
        if ($self->{exonCheck}{isNotFirst}{$startExonID}) {
          print STDERR "WARNING - Do not adjust this start exon (not always first)\n";
        }
        elsif ($extendM or $firstAA ne 'M') {
          # print STDERR "checking upstream region for translation ",$translation->stable_id," because firstAA is $firstAA\n" if $self->{debug};
          my $newCDSStartPos = $self->seekUpstream($extend,$trans,$guide{$translation->stable_id});
          if ($newCDSStartPos) {
            # print STDERR "extending $transStableId CDS upstream $newCDSStartPos\n" if $self->{debug};
            $adjustedStart{$startExonID} = $newCDSStartPos;
            # update the first exon start pos
            my $exonStart = $exons[0]->start;
            my $exonEnd = $exons[0]->end;
            if ($trans->strand == 1) {
              $exonStart = $exons[0]->start($exonStart - $newCDSStartPos - $prepend);
            }
            else {
              $exonEnd = $exons[0]->end($exonEnd + $newCDSStartPos + $prepend);
            }
            $updateFirstExon=1;
          }
        }
      }
      my $endExonID = $trans->end_Exon->dbID;
      if (not exists $adjustedEnd{$endExonID}) {
        if ($self->{exonCheck}{isNotLast}{$endExonID}) {
          print STDERR "WARNING - Do not adjust this end exon (not always last)\n";
        }
        else {
          my $newCDSEndPos = $self->seekDownstream($extend,$trans);
          if ($newCDSEndPos) {
            # print STDERR "extending $transStableId CDS downstream $newCDSEndPos\n" if $self->{debug};
            $adjustedEnd{$endExonID} = $newCDSEndPos;
            # update the last exon end pos
            my $exonStart = $exons[-1]->start;
            my $exonEnd = $exons[-1]->end;
            my $utr_length = $trans->length - $trans->cdna_coding_end;
            if ($trans->strand == 1) {
              $exonEnd = $exons[-1]->end($exonEnd + $newCDSEndPos + $append - $utr_length);
            }
            else {
              $exonStart = $exons[-1]->start($exonStart - $newCDSEndPos - $append + $utr_length);
            }
            $updateLastExon=1;
          }
        }
      }
      if ($updateFirstExon) {
        $changes{gene}{$gene->dbID}=1;
        $changes{trans}{$transID}=1;
        $sth->{update_exon_bounds}->execute($exons[0]->start,$exons[0]->end,$startExonID) unless $self->{nowrite};
        print STDERR "update exon $startExonID ",$exons[0]->start," ",$exons[0]->end,"\n" if $self->{debug};
      }
      if ($updateLastExon and (not $updateFirstExon or @exons > 1)) {
        $changes{gene}{$gene->dbID}=1;
        $changes{trans}{$transID}=1;
        $sth->{update_exon_bounds}->execute($exons[-1]->start,$exons[-1]->end,$endExonID) unless $self->{nowrite};
        print STDERR "update exon $endExonID ",$exons[-1]->start," ",$exons[-1]->end,"\n" if $self->{debug};
      }
    }
    for my $trans (@transcripts) {
      my $transID = $trans->dbID;
      my $transcriptStart = $trans->start;
      my $transcriptEnd = $trans->end;
      next unless $trans->biotype eq 'protein_coding';
      my $translation = $trans->translation;
      $translation or die "huh? I expected to get a translation of transcript " . $trans->stable_id;
      my $translationID = $translation->dbID;
      my $translationStart = $translation->start;
      my $translationEnd = $translation->end;
      my $startExonID = $translation->start_Exon->dbID;
      my $endExonID = $translation->end_Exon->dbID;
      my ($updateTranscript, $updateTranslation);
      if ($adjustedStart{$startExonID}) {
        if ($translationStart > 3) { # some other translation was extended, but this one wasn't
          $translationStart += $adjustedStart{$startExonID};
          $updateTranslation = 1;
        }
        if ($prepend > 0) {
          $translationStart += $prepend;
          $updateTranslation = 1;
        }
        if ($updateTranslation == 1 and $startExonID == $endExonID) {
          $translationEnd += $prepend + $adjustedStart{$startExonID};
        }
        # definitely update start/end of transcript depending on strand
        $updateTranscript=1;
        if ($trans->strand == 1) {
          $transcriptStart = $transcriptStart - $adjustedStart{$startExonID} - $prepend;
        }
        else {
          $transcriptEnd = $transcriptEnd + $adjustedStart{$startExonID} + $prepend;
        }
      }
      if ($adjustedEnd{$endExonID}) {
        $translationEnd += $adjustedEnd{$endExonID};
        $updateTranslation = 1;
    
        # definitely update start/end of transcript depending on strand
        $updateTranscript=1;
        if ($trans->strand == 1) {
          $transcriptEnd = $trans->coding_region_end + $adjustedEnd{$endExonID} + $append;
        }
        else {
          $transcriptStart = $trans->coding_region_start - $adjustedEnd{$endExonID} - $append;
        }
      }
      if ($updateTranscript) {
        $changes{gene}{$gene->dbID}=1;
        $changes{trans}{$transID}=1;
        $sth->{update_transcript_bounds}->execute($transcriptStart,$transcriptEnd,$transID) unless $self->{nowrite};
        print STDERR "update transcript $transID $transcriptStart $transcriptEnd\n" if $self->{debug};
        if ($transcriptStart < $geneStart) {
          $geneStart = $transcriptStart;
          $updateGene=1;
        }
        if ($transcriptEnd > $geneEnd) {
          $geneEnd = $transcriptEnd;
          $updateGene=1;
        }
      }
      if ($updateTranslation) {
        $changes{gene}{$gene->dbID}=1;
        $changes{trans}{$transID}=1;
        $sth->{update_translation_seq}->execute($translationStart,$translationEnd,$translationID) unless $self->{nowrite};
        print STDERR "update translation $translationID $translationStart $translationEnd\n" if $self->{debug};
      }
    }
    if ($updateGene) {
      $sth->{update_gene_bounds}->execute($geneStart,$geneEnd,$geneID) unless $self->{nowrite};
      print STDERR "update gene $geneID $geneStart $geneEnd\n" if $self->{debug};
    }
  }
  return \%changes;
}

sub seekUpstream {
  my $self = shift;
  my $extend = shift;
  my $trans = shift;
  my $estimated_distance = shift;
  $trans->cdna_coding_start <= 3 or return 0;
  my $len = $estimated_distance ? 3*ceil(1.5 * $estimated_distance) : $extend;
  print STDERR "len to check is $len\n" if $self->{debug};
  my $seqToCheck;
  if ($trans->strand == 1) {
    if ($trans->start < $len) {
      $len = $trans->start - ($trans->start % 3) - 1;
      print STDERR "WARNING: + transcript starts near start of seq_region. ",$trans->start," new len to scan is $len\n";
      $len > 2 or return 0;
    }
    my $crs = $trans->coding_region_start;
    $seqToCheck = $trans->slice->subseq($crs-$len, $crs-1, 1);
  }
  else {
    my $cre = $trans->coding_region_end;
    if ($trans->end + $len > $trans->slice->end) {
      $len = $trans->slice->end - $cre;
      $len -= $len % 3;
      print STDERR "WARNING: - transcript starts near end of seq_region. ",$trans->slice->end - $trans->end," new len to scan is $len\n";
      $len > 2 or return 0;
    }
    $seqToCheck = $trans->slice->subseq($cre+1, $cre+$len, -1);
  }
  return findEarliestStart($seqToCheck);
}

sub seekDownstream {
  my $self = shift;
  my $extend = shift;
  my $trans = shift;
  my $utr_length = $trans->length - $trans->cdna_coding_end;
  $utr_length < 3 or return 0;
  my $len = $extend;
  my $seqToCheck;
  if ($trans->strand == 1) {
    my $cre = $trans->coding_region_end;
    if ($cre + $len > $trans->slice->end) {
      $len = $trans->slice->end - $cre;
      $len -= $len % 3;
      print STDERR "WARNING: + transcript ends near end of seq_region. ",$trans->slice->end - $cre," new len to scan is $len\n";
      $len > 2 or return 0;
    }
    $seqToCheck = $trans->slice->subseq($cre-2, $cre+$len, 1);
  }
  else {
    my $crs = $trans->coding_region_start;
    if ($crs < $len) {
      $len = $crs - ($crs % 3) - 1;
      print STDERR "WARNING: - transcript ends near start of seq_region. ",$crs," new len to scan is $len\n";
      $len > 2 or return 0;
    }
    $seqToCheck = $trans->slice->subseq($crs-$len, $crs+2, -1);
  }
  return findFirstStop($seqToCheck);
}

sub findFirstStop {
  my $seq = shift;
  my %stop = (
    TAA => 1,
    TAG => 1,
    TGA => 1
  );
  my $len = length $seq;
  for(my $i=0;$i<$len;$i+=3) {
    my $codon = substr($seq,$i,3);
    return 0 if $codon =~ m/N/;
    return $i if $stop{$codon};
  }
  return 0;
}

sub findEarliestStart {
  my $seq = shift;
  my %codons = (
    TAA => 'stop',
    TAG => 'stop',
    TGA => 'stop',
    ATG => 'start'
  );
  my $len = length $seq;
  my $stop=-1;
  my @starts;
  for(my $i=0;$i<$len;$i+=3) {
    my $codon = substr($seq,$i,3);
    if ($codon =~ m/N/) {
      $stop=$i;
    }
    next unless exists $codons{$codon};
    if ($codons{$codon} eq 'stop') {
      $stop=$i;
    }
    elsif ($codons{$codon} eq 'start') {
      push @starts, $i;
    }
  }
  while (@starts and $starts[0] < $stop) {
    my $skip = shift @starts;
  }
  if (@starts) {
    return $len - $starts[0];
  }
  return 0;
}

sub replaceTranslation {
  my ($self, $params) = @_;
  my %sth = %{$self->{sth}};
  my $minORF = $params->{minORF} || 0;
  my %changes;
  for my $gene (@{$self->{genes}}) {
    next unless $gene->biotype eq 'protein_coding';
    my @transcripts;
    @transcripts = @{$gene->get_all_Transcripts};
    my %updates;
    my $good_translations=0;
    for my $transcript (@transcripts) {
      my $biotype = $transcript->biotype;
      next unless $biotype eq 'protein_coding';
      my $id = $transcript->dbID;
      my $stableid = $transcript->stable_id;
      my $translation = $transcript->translation();
      next unless $translation;

      my $aa_seq;
      eval{$aa_seq= $translation->seq()};
      $@ and die "#ERROR getting translation->seq() for $stableid, $@\n";

      my $bad_translation=0;
      if ($aa_seq) {
        if ($aa_seq !~ /^M/i and $params->{complete}) {
          warn("#WARN aa does not start with M\n");
          $bad_translation=1;
        }
        if ($aa_seq =~ /\*[A-Z]/i) {
          warn("#WARN aa contains internal stop\n");
          $bad_translation=1;
        }
        if ($aa_seq =~ /^x/i or $aa_seq =~ /x$/i) {
          warn("#WARN aa contains leading/trailing x's");
          $bad_translation=1;
        }
      } else {
        warn("#WARN no aa translation for $stableid\n");
        $bad_translation=1;
      }
    
      if (not $bad_translation) {
        $good_translations++;
        next;
      }

      # recalculate translation
      my $mrna = $transcript->spliced_seq();
    
      my @orfs = find_orfs($mrna);
      @orfs or next; 
      my ($orfStart, $orfEnd, $orfLen, $phase) = @{$orfs[0]};
      $changes{gene}{$gene->dbID}=1;
      $changes{trans}{$id}=1;
      if ($orfLen < $minORF) {
        $bad_translation=1;
        # remove translation from db
        print STDERR "delete from translation where transcript_id=",$transcript->dbID,";\n" if $self->{debug};
        $sth{delete_translation}->execute($transcript->dbID) unless $self->{nowrite};
      
        # change transcript from protein-coding to non-coding
        print STDERR "update transcript set biotype='non_coding' where transcript_id=",$transcript->dbID,";\n" if $self->{debug};
        $sth{make_transcript_noncoding}->execute($transcript->dbID) unless $self->{nowrite};
        next;
      }
      $good_translations = 1;
      update_translation($transcript, $orfStart, $orfEnd-3, \%updates, $self->{debug});
    }
    if (not $good_translations) {
      # change gene->biotype to non-coding
      print STDERR "update gene set biotype='non_coding' where gene_id=",$gene->dbID,";\n" if $self->{debug};
      $sth{make_gene_noncoding}->execute($gene->dbID) unless $self->{nowrite};
    }
    $self->doUpdates2(\%updates,\%changes) unless $self->{nowrite};
  }
  return \%changes;
}

sub find_orfs {
    my ( $sequence ) = @_;
    $sequence    = uc $sequence;
    my $codon_table = Bio::Tools::CodonTable->new( -id => 1);
    my $is_start = sub { shift eq 'ATG' };
 
    # stores the begin index of the currently-running ORF in each
    # reading frame
    my @current_orf_start = (-1,-1,-1);
 
    #< stores coordinates of longest observed orf (so far) in each
    #  reading frame
    my @orfs;
 
    # go through each base of the sequence, and each reading frame for each base
    my $seqlen = length $sequence;
    my @start_frame_order;
    for( my $j = 0; $j <= $seqlen-3; $j++ ) {
        my $frame = $j % 3;
 
        my $this_codon = substr( $sequence, $j, 3 );
 
        # if in an orf and this is either a stop codon or the last in-frame codon in the string
        if ( $current_orf_start[$frame] >= 0 ) {
            if ( $codon_table->is_ter_codon( $this_codon ) ||( my $is_last_codon_in_frame = ($j >= $seqlen-5)) ) {
                # record ORF start, end (half-open), length, and frame
                my @this_orf = ( $current_orf_start[$frame], $j+3, undef, $frame );
                my $this_orf_length = $this_orf[2] = ( $this_orf[1] - $this_orf[0] );
                push @orfs, \@this_orf;
                $current_orf_start[$frame] = -1;
            }
        }
        # if this is a start codon
        elsif ( $is_start->($this_codon) ) {
            $current_orf_start[$frame] = $j;
            push @start_frame_order, $frame;
        }
    }
    
    return @orfs ? sort { $b->[2] <=> $a->[2] } @orfs : ();
}

sub update_translation {
  my ($transcript, $orfStart, $orfEnd, $updatesRef, $debug) = @_;
  my $translation = $transcript->translation();
  my @exons = @{$transcript->get_all_Exons};
  my ($seq_start, $seq_end, $start_exon_id, $end_exon_id);
  my $posInTranscript = 0;
  for my $exon (@exons) {
    my ($phase, $end_phase) = (-1, -1);
    if ($orfStart >= $posInTranscript and $orfStart < $posInTranscript + $exon->length) {
      $seq_start = $orfStart - $posInTranscript + 1;
      $start_exon_id = $exon->dbID;
    }
    if ($orfEnd > $posInTranscript and $orfEnd <= $posInTranscript + $exon->length) {
      $seq_end = $orfEnd - $posInTranscript;
      $end_exon_id = $exon->dbID;
    }
    if ($orfStart < $posInTranscript and $posInTranscript < $orfEnd) {
      $phase = ($posInTranscript - $orfStart) % 3;
    }
    $posInTranscript += $exon->length;
    if ($orfStart < $posInTranscript && $posInTranscript <= $orfEnd) {
      $end_phase = ($posInTranscript - $orfStart) % 3;
    }
    # update exon
    print STDERR "update exon set phase=$phase, end_phase=$end_phase where exon_id=",$exon->dbID,";\n" if $debug;
    # $sth{update_exon}->execute($phase, $end_phase, $exon->dbID) unless $nowrite;

    $updatesRef->{$exon->dbID}{$transcript->dbID} = {
      exonID => $exon->dbID, # this could change if a new exon needs to be created
      phase => $phase,
      end_phase => $end_phase
    };
  }
  # update translation
  print STDERR "update translation set start_exon_id=$start_exon_id, seq_start=$seq_start, end_exon_id=$end_exon_id, seq_end=$seq_end where translation_id=",$translation->dbID,";\n" if $debug;
  # $sth{update_translation}->execute($start_exon_id, $seq_start, $end_exon_id, $seq_end, $translation->dbID) unless $nowrite;
  my $startUpdate = $updatesRef->{$start_exon_id}{$transcript->dbID};
  $startUpdate->{translationID} = $translation->dbID;
  $startUpdate->{seq_start} = $seq_start;
  my $endUpdate = $updatesRef->{$end_exon_id}{$transcript->dbID};
  $endUpdate->{translationID} = $translation->dbID;
  $endUpdate->{seq_end} = $seq_end;
}

sub doUpdates {
  my ($self, $updatesRef) = @_;
  print STDERR "doUpdates ",Dumper($updatesRef) if $self->{debug};
  my %updates = %$updatesRef;
  for my $eid (keys %updates) {
    # check if updates are consistent w.r.t. exon phase
    # if not, create a new exon for the updated translation
    my %phases;
    my $orig_phase;
    for my $tid (keys %{$self->{exon_transcript_phase}{$eid}}) {
      $orig_phase = $self->{exon_transcript_phase}{$eid}{$tid};
      if (exists $updates{$eid}{$tid}) {
        $phases{$updates{$eid}{$tid}{phase}}{$tid} = 1;
      }
      else {
        $phases{$orig_phase}{$tid} = 1;
      }
    }
    if (keys %phases > 1) {
      # create a new exon for each new phase
      $self->{sth}{select_exon}->execute($eid);
      my ($eid2,$sr,$st,$en,$str,$ph,$ep,$cur,$con,$sid,$vers,@dates) = @{$self->{sth}{select_exon}->fetchrow_arrayref};
      for my $phase (keys %phases) {
        $phase = $phase + 0; # hash keys are strings, but phase is an integer
        next if $phase == $orig_phase;
        print STDERR "creating new exon with phase $phase\n" if $self->{debug};
        $vers++;
        $self->{sth}{insert_exon}->execute(undef,$sr,$st,$en,$str,$phase,$ep,$cur,$con,$sid,$vers);
        my $newExonId = $self->{dbh}->last_insert_id(undef, undef, undef, undef);
        print STDERR "got new exon id $newExonId\n" if $self->{debug};
        # update exon_transcript table and startExon of translation (update happens after this if block)
        for my $tid (keys %{$phases{$phase}}) {
          print STDERR "updating exon_transcript $newExonId,$eid,$tid\n" if $self->{debug};
          $self->{sth}{update_exon_transcript}->execute($newExonId,$eid,$tid);
          $self->{exon_transcript_phase}{$newExonId}{$tid} = $phase;
          $updates{$eid}{$tid}{endExon} = $newExonId if ($updates{$eid}{$tid}{startExon} == $updates{$eid}{$tid}{endExon});
          $updates{$eid}{$tid}{startExon} = $newExonId;
        }
      }
    }
    else {
      my ($phase) = keys %phases;
      $phase = $phase + 0;
      if ($phase != $orig_phase) {
        print STDERR "updating exon phase $eid $phase\n" if $self->{debug};
        $self->{sth}{update_exon_phase}->execute($phase, $eid);
        for my $tid (keys %{$self->{exon_transcript_phase}{$eid}}) {
          $self->{exon_transcript_phase}{$eid}{$tid} = $phase;
        }
      }
    }
    # update the translations
    for my $tid (keys %{$updates{$eid}}) {
      my $u = $updates{$eid}{$tid};
      print STDERR "updating translation ",join(',',$u->{startExon},$u->{seqStart},$u->{endExon},$u->{translationID}),"\n" if $self->{debug};
      $self->{sth}{update_translation}->execute($u->{startExon}, $u->{seqStart}, $u->{endExon}, $u->{translationID});
    }
  }
}

sub doUpdates2 {
  my ($self, $updatesRef) = @_;
  my %updates = %$updatesRef;
  for my $eid (keys %updates) {
    # check if updates are consistent w.r.t. exon phase
    # if not, create a new exon for the updated translation
    my %phases;
    my $orig_phase;
    for my $tid (keys %{$self->{exon_transcript_phase}{$eid}}) {
      $orig_phase = $self->{exon_transcript_phase}{$eid}{$tid};
      if (exists $updates{$eid}{$tid}) {
        $phases{$updates{$eid}{$tid}{phase}}{$tid} = 1;
      }
      else {
        $phases{$orig_phase}{$tid} = 1;
      }
    }
    if (keys %phases > 1) {
      # create a new exon for each new phase
      $self->{sth}{select_exon}->execute($eid);
      my ($eid2,$sr,$st,$en,$str,$ph,$ep,$cur,$con,$sid,$vers,@dates) = @{$self->{sth}{select_exon}->fetchrow_arrayref};
      for my $phase (keys %phases) {
        $phase = $phase + 0; # hash keys are strings, but phase is an integer
        next if $phase == $orig_phase;
        print STDERR "creating new exon with phase $phase\n" if $self->{debug};
        $vers++;
        $self->{sth}{insert_exon}->execute(undef,$sr,$st,$en,$str,$phase,$ep,$cur,$con,$sid,$vers);
        my $newExonId = $self->{dbh}->last_insert_id(undef, undef, undef, undef);
        print STDERR "got new exon id $newExonId\n" if $self->{debug};
        # update exon_transcript table and startExon of translation (update happens after this if block)
        for my $tid (keys %{$phases{$phase}}) {
          print STDERR "updating exon_transcript $newExonId,$eid,$tid\n" if $self->{debug};
          $self->{sth}{update_exon_transcript}->execute($newExonId,$eid,$tid);
          $self->{exon_transcript_phase}{$newExonId}{$tid} = $phase;
          $updates{$eid}{$tid}{exonID} = $newExonId;
        }
      }
    }
    else {
      my ($phase) = keys %phases;
      $phase = $phase + 0;
      if ($phase != $orig_phase) {
        print STDERR "updating exon phase $eid $phase\n" if $self->{debug};
        $self->{sth}{update_exon_phase}->execute($phase, $eid);
        for my $tid (keys %{$self->{exon_transcript_phase}{$eid}}) {
          $self->{exon_transcript_phase}{$eid}{$tid} = $phase;
        }
      }
    }
    # update the translations
    for my $tid (keys %{$updates{$eid}}) {
      my $u = $updates{$eid}{$tid};
      if ($u->{seq_start}) {
        print STDERR "updating translation start\n" if $self->{debug};
        $self->{sth}{update_translation_start}->execute($u->{exonID}, $u->{seq_start}, $u->{translationID});
      }
      if ($u->{seq_end}) {
        print STDERR "updating translation end\n" if $self->{debug};
        $self->{sth}{update_translation_end}->execute($u->{exonID}, $u->{seq_end}, $u->{translationID});
      }
    }
  }
}

1;

__END__

