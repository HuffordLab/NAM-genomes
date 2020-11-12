#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Data::Dumper;

my ($gene_models_file,$interproscan_output, $MIN_TPM, $MIN_EVI_COVERAGE, $MIN_REF_COVERAGE, $MAX_AED, @stringtie_files) = @ARGV;
($MAX_AED >= 0 and $MAX_AED <= 1) or die "usage: $0 gff ipr minTPM min_evi_cov min_ref_cov max_aed stringtie_files\n";
print STDERR "parsing gene models\n";
my $models = parseGFF($gene_models_file);
print STDERR "loading interproscan results\n";
my $domains = parse_interproscan($interproscan_output);
my %samples;
print STDERR "loading stringtie files\n";
for my $gtf (@stringtie_files) {
  if ($gtf =~ m/(\S+)\.gtf$/) {
    print STDERR "$1\n";
    $samples{$1} = parseStringtie($gtf);
  }
}
print STDERR "the magic begins\n";

for my $chr (sort keys %$models) {
  for my $strand (sort keys %{$models->{$chr}}) {
    my @mod_genes = sort {$a->{start} <=> $b->{start}} @{$models->{$chr}{$strand}};
    for my $sample (keys %samples) {
      my $evidence = $samples{$sample};
      my @evi_genes = ();
      if ($evidence->{$chr} and $evidence->{$chr}{$strand}) {
        @evi_genes = sort {$a->{start} <=> $b->{start}} @{$evidence->{$chr}{$strand}};
      }
      score_overlaps(\@mod_genes, \@evi_genes, $sample);
    }
    for my $gene (@mod_genes) {
      my ($winners,$score) = rank_transcripts($gene, $domains);
      for (my $i=0;$i<@$winners;$i++) {
        my $tid = $winners->[$i];
        my $tr = $gene->{transcripts}{$tid};
        my $domain_len = $domains->{$tid} || 0;
        print join("\t",$gene->{gene_id}, $tid, $i+1, $score->{$tid}, $tr->{len}, $tr->{cds}, $domain_len, map {$tr->{AED}{$_}} sort keys %{$tr->{AED}}),"\n";
      }
    }
  }
}
sub rank_transcripts {
  my ($gene, $domains) = @_;
  my @candidates = keys %{$gene->{transcripts}};
  # gather rankings for each sample (by aed) and length metric
  my %aeds;
  for my $tid (@candidates) {
    my $tr = $gene->{transcripts}{$tid};
    $tr->{domain} = $domains->{$tid} || 0;
    for my $sample (keys %{$tr->{AED}}) {
      if ($tr->{AED}{$sample} <= $MAX_AED) {
        if (not exists $aeds{$sample}) {
          $aeds{$sample} = [];
        }
        push @{$aeds{$sample}}, {
          tid => $tid,
          aed => $tr->{AED}{$sample}
        }
      }
    }
  }
  # rank transcripts for each sample
  my %rank;
  my %score; # sum of weighted votes for each transcript where weight = 1 - AED
  # initialize score hash with length metrics
  # longest domain is 100, weight for tid with domain 50 is 0.5
  # if longest domain is 0, initialize to 0
  my %max = (domain => 0, cds => 0, len => 0);
  my %metricWeight = (domain => 3, cds => 2, len => 1);
  for my $tid (@candidates) {
    $score{$tid} = 0;
    my $tr = $gene->{transcripts}{$tid};
    for my $metric (keys %max) {
      $max{$metric} = $tr->{$metric} if ($tr->{$metric} > $max{$metric});
    }
  }
  for my $tid (@candidates) {
    my $tr = $gene->{transcripts}{$tid};
    for my $metric (keys %max) {
      if ($max{$metric} > 0) {
        $score{$tid} += $metricWeight{$metric} * $tr->{$metric}/$max{$metric};
      }
    }
  }
  for my $sample (keys %aeds) {
    my @sorted = sort {$a->{aed} <=> $b->{aed}} grep {$_->{aed} <= $MAX_AED} @{$aeds{$sample}};
    my $prevAED = -1;
    my $rank_i=0;
    for my $ta (@sorted) {
      $score{$ta->{tid}} += 1 - $ta->{aed};
      if ($ta->{aed} > $prevAED) {
        $prevAED = $ta->{aed};
        $rank_i++;
      }
      $rank{$sample}{$ta->{tid}} = $rank_i;
    }
  }
  # rank by length
  my @sorted_by_length = sort {
    $gene->{transcripts}{$b}{domain} <=> $gene->{transcripts}{$a}{domain} or
    $gene->{transcripts}{$b}{cds} <=> $gene->{transcripts}{$a}{cds} or
    $gene->{transcripts}{$b}{len} <=> $gene->{transcripts}{$a}{len}
  } @candidates;
  my $length_rank=0;
  my $prev = {domain => 999999};
  for my $tid (@sorted_by_length) {
    my $tr = $gene->{transcripts}{$tid};
    if (
      $tr->{domain} < $prev->{domain} or
      $tr->{cds} < $prev->{cds} or
      $tr->{len} < $prev->{len}
    ) {
      $prev = $tr;
      $length_rank++;
    }
    $rank{_length}{$tid} = $length_rank;
  }
  # print STDERR Dumper(\%rank);
  my @winners;
  my $nTranscripts = @sorted_by_length;
  while (@candidates > 1) {
    my %votes;
    for my $voter (keys %rank) {
      my @ranked = sort {$rank{$voter}{$a} <=> $rank{$voter}{$b}} keys %{$rank{$voter}};
      @ranked or next;
      my $prevRank = $rank{$voter}{$ranked[0]}-1;
      my $pri=$prevRank-1;
      for my $tid (@ranked) {
        # print STDERR "$tid $prevRank $pri\n";
        if ($rank{$voter}{$tid} > $prevRank) {
          $pri++;
          $prevRank = $rank{$voter}{$tid};
        }
        # print STDERR "vote for $tid in $pri\n";
        $votes{$pri}{$tid}++;
      }
    }
    # print STDERR "got all votes\n";
    # print STDERR Dumper(\%votes);
    my $round = 0;
    my %inTheRunning;
    for my $tid (@candidates) {
      $inTheRunning{$tid}=1;
    }
    my $leader;
    my $maxRounds = scalar keys %votes;
    while (keys %inTheRunning > 1 and $round < $maxRounds) {
      # print STDERR "round $round of $maxRounds inTheRunning ",join(',',keys %inTheRunning),"\n";
      my @sort_by_votes = sort {$votes{$round}{$b} <=> $votes{$round}{$a}} grep {exists $votes{$round}{$_}} keys %inTheRunning;
      if (@sort_by_votes) {
        %inTheRunning = ();
        $leader = shift @sort_by_votes;
        $inTheRunning{$leader} = 1;
        # print STDERR "leader $leader has $votes{$round}{$leader} votes\n";
        for my $tid (@sort_by_votes) {
          # print STDERR "$tid $round $votes{$round}{$tid}\n";
          last if ($votes{$round}{$tid} < $votes{$round}{$leader});
          $inTheRunning{$tid}=1;
        }
      }
      $round++;
    }
    my @results = sort keys %inTheRunning;
    # print STDERR "After $round rounds of voting: [@results]\n";
    
    my $winner = shift @results;
    push @winners, $winner;
    
    # remove winner from consideration in next round
    for my $voter (keys %rank) {
      delete $rank{$voter}{$winner} if exists $rank{$voter}{$winner};
    }
    @candidates = keys %{$rank{_length}};
  }
  # output candidate
  push @winners, shift @candidates;
  return (\@winners,\%score);
}

sub parseGFF {
  my $filename = shift;
  my %span;
  my %transcripts; # key is gene_id value is hash ref with key transcript_id
  my %genes; # keys are {chr}{strand} values are gene_ids
  my %g2t;
  open (my $fh, "<", $filename);
  while (<$fh>) {
    next if (/^#/);
    chomp;
    my ($chr,$method,$type,$start,$end,$score,$strand,$idk,$attr_string) = split /\t/, $_;
    $chr =~ s/^chr//;
    my %attr;
    for my $kv (split /;/, $attr_string) {
      my ($k, $v) = split /=/, $kv;
      if ($k eq '_QI') {
        $attr{$k} = split /\|/, $v;
      }
      elsif ($k eq 'Dbxref' or $k eq 'Ontology_term') {
        $attr{$k} = split /,/, $v;
      }
      else {
        $v =~ s/^gene://;
        $v =~ s/^transcript://;
        $attr{$k} = $v;
      }
    }
    if ($type eq 'gene') {
      $span{$attr{ID}} = [$start,$end];
      if (not exists $genes{$chr}{$strand}) {
        $genes{$chr}{$strand} = [];
      }
      push @{$genes{$chr}{$strand}}, $attr{ID};
    }
    if ($type eq 'mRNA') {
      my $transcript_id = $attr{ID};
      my $gene_id = $attr{Parent};
      if (not exists $transcripts{$transcript_id}) {
        $transcripts{$transcript_id} = {exons=>[],len=>0, cds=>0};
      }
      $transcripts{$transcript_id}{attr} = \%attr;
      $g2t{$gene_id}{$transcript_id} = $transcripts{$transcript_id};
    }
    if ($type eq 'exon') {
      for my $transcript_id (split /,/, $attr{Parent}) {
        if (not exists $transcripts{$transcript_id}) {
          $transcripts{$transcript_id} = {exons=>[],len=>0, cds=>0};
        }
        push @{$transcripts{$transcript_id}{exons}}, [$start,$end];
        $transcripts{$transcript_id}{len} += $end - $start + 1;
      }
    }
    if ($type eq 'CDS') {
      my $transcript_id = $attr{Parent};
      $transcripts{$transcript_id}{cds} += $end - $start + 1;
    }
  }
  close $fh;
  # sort genes by span
  for my $chr (keys %genes) {
    for my $strand (keys %{$genes{$chr}}) {
      for (my $i = 0; $i<@{$genes{$chr}{$strand}}; $i++) {
        my $gene_id = $genes{$chr}{$strand}[$i];
        $genes{$chr}{$strand}[$i] = {
          gene_id => $gene_id,
          start   => $span{$gene_id}[0],
          end     => $span{$gene_id}[1],
          transcripts => $g2t{$gene_id}
        };
      }
    }
  }
  return \%genes;
}

sub parseStringtie {
  my $filename = shift;
  my %span; # key is gene_id
  my %transcripts; # key is gene_id value is hash ref with key transcript_id
  my %genes; # keys are {chr}{strand} values are gene_ids
  open (my $fh, "<", $filename);
  while (<$fh>) {
    next if (/^#/);
    chomp;
    my ($chr,$method,$type,$start,$end,$score,$strand,$idk,$attribute_string) = split /\t/, $_;
    $chr =~ s/^chr//;
    my ($gene_id, $transcript_id) = $attribute_string =~ m/gene_id\s+"(\S+)";\s+transcript_id\s+"(\S+)";/;
    my ($tpm) = $attribute_string =~ m/;\sTPM\s"(\S+)";/;
    if ($type eq 'transcript') {
      $transcripts{$gene_id}{$transcript_id} = {exons=>[],len=>0};
      if ($tpm) {
        $transcripts{$gene_id}{$transcript_id}{tpm} = $tpm;
      }
      if (not exists $span{$gene_id}) {
        $span{$gene_id} = [$start,$end];
        if (not exists $genes{$chr}{$strand}) {
          $genes{$chr}{$strand} = [];
        }
        push @{$genes{$chr}{$strand}}, $gene_id;
      }
    }
    if ($type eq 'exon') {
      if (@{$transcripts{$gene_id}{$transcript_id}{exons}}
      and $start < $transcripts{$gene_id}{$transcript_id}{exons}[0][0]) {
        unshift @{$transcripts{$gene_id}{$transcript_id}{exons}}, [$start,$end];
      }
      else {
        push @{$transcripts{$gene_id}{$transcript_id}{exons}}, [$start,$end];
      }
      $transcripts{$gene_id}{$transcript_id}{len} += $end - $start + 1;
      if ($start < $span{$gene_id}[0]) {
        $span{$gene_id}[0] = $start;
      }
      if ($end > $span{$gene_id}[1]) {
        $span{$gene_id}[1] = $end;
      }
    }
  }
  close $fh;

  # sort genes by span
  for my $chr (keys %genes) {
    for my $strand (keys %{$genes{$chr}}) {
      for (my $i = 0; $i<@{$genes{$chr}{$strand}}; $i++) {
        my $gene_id = $genes{$chr}{$strand}[$i];
        $genes{$chr}{$strand}[$i] = {
          gene_id => $gene_id,
          start   => $span{$gene_id}[0],
          end     => $span{$gene_id}[1],
          transcripts => $transcripts{$gene_id}
        };
      }
    }
  }
  return \%genes;
}

sub score_overlaps {
  my ($ref, $evi, $sample) = @_;
  # ref and evi are references to sorted arrays of genes
  # each gene has the fields: gene_id, start, end, and transcripts
  # transcripts is a hashref with transcript_id keys and fields: length and exons
  # exons is a [start,end]
  
  # iterate over the ref genes and the evi genes
  # when a pair of genes overlaps by at least threshold percent
  # calculate a scoring matrix for the transcripts of the ref gene vs the transcripts of the evi gene
  my %minDist;
  my %maxTPMDist;
  my $evi_offset = 0;
  my $evi_length = @$evi;
  for (my $ref_offset = 0; $ref_offset < @$ref; $ref_offset++) {
    # initialize distances for each transcript
    my @transcripts = sort keys %{$ref->[$ref_offset]{transcripts}};
    for my $ref_tid (@transcripts) {
      $ref->[$ref_offset]{transcripts}{$ref_tid}{AED}{$sample} = [];
    }
    while ($evi_offset < $evi_length - 1
    and $evi->[$evi_offset]{end} < $ref->[$ref_offset]{start}) {
      $evi_offset++;
    }
    my $evi_first = $evi_offset;
    while ($evi_offset < $evi_length
    and $evi->[$evi_offset]{end} > $ref->[$ref_offset]{start}
    and $evi->[$evi_offset]{start} < $ref->[$ref_offset]{end}) {
      # genes overlap
      my $aed = score_transcripts2($ref->[$ref_offset], $evi->[$evi_offset]);
      for my $ref_tid (@transcripts) {
        if (exists $aed->{$ref_tid}) {
          push @{$ref->[$ref_offset]{transcripts}{$ref_tid}{AED}{$sample}}, $aed->{$ref_tid};
        }
      }
      $evi_offset++;
    }
    $evi_offset = $evi_first;
    # average the AEDs we just calculated
    for my $ref_tid (@transcripts) {
      my @aeds = @{$ref->[$ref_offset]{transcripts}{$ref_tid}{AED}{$sample}};
      if (@aeds) {
        my $sum=0;
        for my $aed (@aeds) {
          $sum += $aed;
        }
        $ref->[$ref_offset]{transcripts}{$ref_tid}{AED}{$sample} =  $sum / scalar @aeds;
      }
      else {
        $ref->[$ref_offset]{transcripts}{$ref_tid}{AED}{$sample} = 1;
      }
    }
  }
}

sub score_transcripts2 {
  my ($ref, $evi, $maxTPMDist) = @_;
  # find the dominant transcript, measure the overlap, possibly trim the ref transcript
  my @sort_by_tpm = sort {$evi->{transcripts}{$b}{tpm} <=> $evi->{transcripts}{$a}{tpm}} grep {
    $evi->{transcripts}{$_}{tpm} >= $MIN_TPM
  } keys %{$evi->{transcripts}};
  return {} unless @sort_by_tpm;
  my $maxTPM_tid = $sort_by_tpm[0];
  my $evi_transcript = $evi->{transcripts}{$maxTPM_tid};
  my $evi_start = $evi_transcript->{exons}[0][0];
  my $evi_end = $evi_transcript->{exons}[-1][1];
  my %myAED;
  for my $ref_tid (sort keys %{$ref->{transcripts}}) {
    my $ref_transcript = $ref->{transcripts}{$ref_tid};
    my $ref_start = $ref_transcript->{exons}[0][0];
    my $ref_end = $ref_transcript->{exons}[-1][1];
    if (@{$ref_transcript->{exons}[0]} == 0) {
      print STDERR "no exons in $ref_tid\n";
      next;
    }
    my $exon_overlap = measure_overlap($ref_transcript->{exons}, $evi_transcript->{exons});
    if ($exon_overlap/$evi_transcript->{len} > $MIN_EVI_COVERAGE
    && $exon_overlap/$ref_transcript->{len} > $MIN_REF_COVERAGE) {
      my $adjusted_ref_len = $ref_transcript->{len};
      my $adjusted_evi_len = $evi_transcript->{len};
      if ($ref_start < $evi_start) {
        $adjusted_ref_len -= measure_overlap($ref_transcript->{exons}, [[$ref_start,$evi_start-1]]);
      }
      if ($ref_start > $evi_start) {
        $adjusted_evi_len -= measure_overlap($evi_transcript->{exons}, [[$evi_start-1,$ref_start-1]])
      }
      if ($ref_end > $evi_end) {
        $adjusted_ref_len -= measure_overlap($ref_transcript->{exons}, [[$evi_end+1, $ref_end]]);
      }
      if ($ref_end < $evi_end) {
        $adjusted_evi_len -= measure_overlap($evi_transcript->{exons}, [[$ref_end+1, $evi_end]]);
      }
      if ($adjusted_ref_len <= 0) {
        die "$ref_tid adjusted_ref_len $adjusted_ref_len <= 0\n";
      }
      if ($adjusted_evi_len <= 0) {
        die "$ref_tid adjusted_evi_len $adjusted_evi_len <= 0\n";
      }
      my $D = 1 - ($exon_overlap/$adjusted_ref_len + $exon_overlap/$adjusted_evi_len)/2;
      $myAED{$ref_tid} = $D;
    }
  }
  return \%myAED;
}

sub measure_overlap {
  my ($t1, $t2) = @_;
  
  my %edges;
  for my $exon (@$t1, @$t2) {
    my ($start, $end) = @$exon;
    $edges{$start} = exists $edges{$start} ? $edges{$start}+1 : 1;
    $edges{$end} = exists $edges{$end} ? $edges{$end}-1 : -1;
  }
  my $depth = 0;
  my $since = 0;
  my $overlap = 0;
  for my $pos (sort {$a <=> $b} keys %edges) {
    # print STDERR "MO $pos $edges{$pos}\n";
    $depth += $edges{$pos};
    if ($depth == 2) {
      $since = $pos;
    }
    elsif ($since > 0) {
      $overlap += $pos - $since + 1;
      $since = 0;
    }
  }
  return $overlap;
}


sub parse_interproscan {
  my $ipr_file = shift;
  my %domains;
  open (my $fh, "<", $ipr_file);
  while (<$fh>) {
    chomp;
    my ($tid,$x,$len,$db,$id,$desc,$start,$end,$score,$t,$date,$ipr,$desc2,$go) = split /\t/, $_;
    if ($ipr) {
      $tid =~ s/_P/_T/;
      $start--;
      exists $domains{$tid}{$start} or $domains{$tid}{$start} = 0;
      exists $domains{$tid}{$end} or $domains{$tid}{$end} = 0;
      $domains{$tid}{$start}++;
      $domains{$tid}{$end}--;
    }
  }
  close $fh;
  for my $tid (keys %domains) {
    my $depth=0;
    my $len=0;
    my $start=0;
    for my $pos (sort {$a <=> $b} keys %{$domains{$tid}}) {
      if ($depth == 0) {
        # beginning of span
        $start = $pos;
      }
      $depth += $domains{$tid}{$pos};
      if ($depth == 0) {
        # end of span
        $len += $pos - $start;
      }
    }
    $domains{$tid} = $len;
  }
  return \%domains;
}
