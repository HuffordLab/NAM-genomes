#!/usr/bin/perl
#
#  FASTA Splitter  -  a script for partitioning a FASTA file into pieces
#
#  Version 0.2.0 (February 14, 2014)
#
#  Copyright (c) 2012-2014 Kirill Kryukov
#
#  This software is provided 'as-is', without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.
#

use Getopt::Long;
use strict;
$| = 1;

my $start_time = time;

my ($opt_n_parts,$opt_part_size,$opt_measure,$opt_line_len,$opt_eol,$opt_version,$opt_help);
GetOptions("n-parts=i"     => \$opt_n_parts,
           "part-size=i"   => \$opt_part_size,
           "measure=s"     => \$opt_measure,
           "line-length=i" => \$opt_line_len,
           "eol=s"         => \$opt_eol,
           "version"       => \$opt_version,
           "help"          => \$opt_help)
or die "Can't parse command line arguments\n";

if ($opt_version) { show_version(); }
if ($opt_help) { show_help(); }

if (!defined($opt_n_parts) and !defined($opt_part_size))
{
    if ($opt_version || $opt_help) { exit(0); }
    else { die "Splitting method is not specified\nUse -h for help\n"; }
}
if (!@ARGV) { die "File for splitting is not specified\n"; }

if (defined($opt_n_parts) and $opt_n_parts <= 0) { die "Non-positive number of parts\n"; }
if (defined($opt_part_size) and $opt_part_size <= 0) { die "Non-positive part size\n"; }
if (defined($opt_measure) and $opt_measure ne 'all' and $opt_measure ne 'seq' and $opt_measure ne 'count') { die "Unknown value of --measure option\n"; }
if (defined($opt_eol) and $opt_eol ne 'dos' and $opt_eol ne 'mac' and $opt_eol ne 'unix') { die "Unknown value of --eol option\n"; }

my $n_parts = defined($opt_n_parts) ? $opt_n_parts : 0;
my $part_size = defined($opt_part_size) ? $opt_part_size : 0;
my $line_len = (defined($opt_line_len) and $opt_line_len >= 0) ? $opt_line_len : 60;
my $eol = defined($opt_eol) ? (($opt_eol eq 'dos') ? "\x0D\x0A" : ($opt_eol eq 'mac') ? "\x0D" : "\x0A") : "\x0A";
my $eol_len = length($eol);
my $measure = defined($opt_measure) ? (($opt_measure eq 'count') ? 0 : ($opt_measure eq 'seq') ? 1 : 2) : 2;
my @part_start = ();
my ($base,$ext,$num_len,$total_size);
my ($OUT,$name,$data,$written_total,$written_this_part,$part_end,$part);

foreach my $infile (@ARGV) { split_file($infile); }

my $end_time = time;
my $elapsed_time = $end_time - $start_time;
print "All done, $elapsed_time second", (($elapsed_time==1)?'':'s'), " elapsed\n";

sub split_file
{
    my ($infile) = @_;
    if (!-e $infile or !-f $infile) { print "Can't find file \"$infile\"\n"; return; }
    print $infile;

    ($base,$ext) = ($infile,'');
    if ($base =~ /[\/\\]([^\/\\]+)$/) { $base = $1; }
    if ($base =~ /^(.+?)(\.(fasta|faa|fna|fa))$/i) { ($base,$ext) = ($1,$2); }

    @part_start = ();
    my ($n_seq,$total_seq_len,$n_parts_found) = (0,0,0);

    if ($part_size)
    {
        ($n_seq,$total_seq_len,$total_size,$n_parts_found) = get_file_size_and_part_boundaries($infile);
        if (!$n_parts) { print ": $n_seq sequences, $total_seq_len bp"; }
        print ' => ', ($n_parts ? 'extracting' : 'dividing into'), ' ', $n_parts_found, ' part', ($n_parts_found > 1 ? 's' : ''),
              " of <= $part_size ", ($measure ? (($measure > 1) ? 'bytes' : 'bp') : 'sequences'), "\n";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        binmode $IN;
        $num_len = length($n_parts_found);
        $OUT = undef;
        my ($out_file,$part,$si,$buffer) = (undef,0,-1,'');
        while (<$IN>)
        {
            $_ =~ s/[\x0D\x0A]+$//;
            if (substr($_,0,1) eq '>')
            {
                if ($OUT)
                {
                    if ($line_len == 0) { if ($si > 1) { print $OUT $eol; } }
                    elsif ($buffer ne '') { print $OUT $buffer, $eol; $buffer = ''; }
                }
                $si++;
                if ($si >= $part_start[$part+1])
                {
                    if ($OUT) { close $OUT; }
                    $part++;
                    if ($part > $n_parts_found) { last; }
                    $out_file = sprintf("%s.part-%0*d%s",$base,$num_len,$part,$ext);
                    open($OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
                    binmode $OUT;
                }
                print $OUT $_, $eol;
                next;
            }
            if (!$line_len) { print $OUT $_; }
            else
            {
                $buffer .= $_;
                while (length($buffer) >= $line_len) { print $OUT substr($buffer,0,$line_len,''), $eol; }
            }
        }
        close $IN;
        if ($OUT)
        {
            if (!$line_len) { if ($si > 1) { print $OUT $eol; } }
            elsif ($buffer ne '') { print $OUT $buffer, $eol; $buffer = ''; }
            close $OUT;
        }
    }
    else
    {
        ($n_seq,$total_seq_len,$total_size) = get_file_size($infile);
        print ": $n_seq sequences, $total_seq_len bp => dividing into $n_parts part", ($n_parts > 1 ? 's' : ''), " ";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        binmode $IN;
        $num_len = length($n_parts);
        ($OUT,$name,$data,$written_total,$written_this_part,$part_end,$part) = (undef,undef,'',0,0,int($total_size / $n_parts),1);
        while(<$IN>)
        {
            $_ =~ s/[\x0D\x0A]+$//;
            if (substr($_,0,1) eq '>')
            {
                if (defined $name) { dump_seq(); }
                $name = $_; $data = ''; next;
            }
            $data .= $_;
        }
        if (defined $name) { dump_seq(); }
        close $IN;
        if ($OUT) { close $OUT; }
        print " OK\n";
    }
}

sub dump_seq
{
    my $slen = length($data);
    my $seq_size = ($measure == 0) ? 1 : ($measure == 1) ? $slen : $slen + length($name) + $eol_len*(1 + ($line_len ? int(($slen+$line_len-1)/$line_len) : 1));
    my $new_written_total = $written_total + $seq_size;
    if ( !$OUT or
         ($written_this_part and ($new_written_total > $part_end) and ($new_written_total - $part_end > $part_end - $written_total)) )
    {
        if ($OUT) { close $OUT; $part++; }
        my $part_file = $base;
        if ($part_file !~ /\.part-\d+$/) { $part_file .= '.part'; }
        $part_file .= sprintf("-%0*d%s",$num_len,$part,$ext);
        open($OUT,'>',$part_file) or die "Error: Can't create file \"$part_file\"\n";
        binmode $OUT;
        $part_end = int($total_size / $n_parts * $part) + 1;
        $written_this_part = 0;
        print ".";
    }
    print $OUT $name, $eol;
    if ($line_len) { for (my $s=0; $s<$slen; $s+=$line_len) { print $OUT substr($data,$s,$line_len), $eol; } }
    else { print $OUT $data, $eol; }
    $written_this_part += $seq_size;
    $written_total += $seq_size;
}

sub get_file_size_and_part_boundaries
{
    my ($file) = @_;
    open(my $IN,'<',$file) or die "Error: Can't open file \"$file\"\n";
    binmode $IN;
    my ($nseq,$total_seq_length,$total_size,$n_parts_found,$this_part_size,$nlen,$slen,$stop) = (0,0,0,1,0,0,0,0);
    $part_start[1] = 0;
    while (<$IN>)
    {
        $_ =~ s/[\x0D\x0A]+$//;
        my $len = length($_);
        if (substr($_,0,1) eq '>')
        {
            if ($nlen)
            {
                my $seq_size = seq_size($nlen,$slen);
                if ($part_size and $this_part_size and ($this_part_size + $seq_size > $part_size))
                {
                    if ($n_parts and $n_parts_found == $n_parts) { $stop = 1; last; }
                    else { $this_part_size = $seq_size; $n_parts_found++; $part_start[$n_parts_found] = $nseq; }
                }
                else { $this_part_size += $seq_size; }
                $nseq++; $total_seq_length += $slen; $total_size += $seq_size; 
            }
            ($nlen,$slen) = ($len,0); next;
        }
        if ($nlen) { $slen += $len; }
    }
    if ($nlen and !$stop)
    {
        my $seq_size = seq_size($nlen,$slen);
        if ($part_size and $this_part_size and ($this_part_size + $seq_size > $part_size))
        {
            if ($n_parts and $n_parts_found == $n_parts) { $stop = 1; }
            else { $this_part_size = $seq_size; $n_parts_found++; $part_start[$n_parts_found] = $nseq; }
        }
        if (!$stop) { $nseq++; $total_seq_length += $slen; $total_size += $seq_size; }
    }
    close $IN;
    $part_start[$n_parts_found+1] = $nseq;
    return ($nseq,$total_seq_length,$total_size,$n_parts_found);
}

sub get_file_size
{
    my ($file) = @_;
    open(my $IN,'<',$file) or die "Error: Can't open file \"$file\"\n";
    binmode $IN;
    my ($nseq,$total_seq_length,$total_size,$nlen,$slen) = (0,0,0,0,0);
    while (<$IN>)
    {
        $_ =~ s/[\x0D\x0A]+$//;
        my $len = length($_);
        if (substr($_,0,1) eq '>')
        { 
            if ($nlen) { $nseq++; $total_seq_length += $slen; $total_size += seq_size($nlen,$slen); }
            ($nlen,$slen) = ($len,0); next;
        }
        if ($nlen) { $slen += $len; }
    }
    if ($nlen) { $nseq++; $total_seq_length += $slen; $total_size += seq_size($nlen,$slen); }
    close $IN;
    return ($nseq,$total_seq_length,$total_size);
}

sub seq_size
{
    my ($nlen,$slen) = @_;
    return ($measure == 0) ? 1 :
           ($measure == 1) ? $slen :
           $slen + $nlen + $eol_len*(1 + ($line_len ? int(($slen+$line_len-1)/$line_len) : 1));
}

sub show_version
{
    print q{FASTA Splitter 0.2.0 (February 14, 2014)
Copyright (c) 2012-2014 Kirill Kryukov
};
}

sub show_help()
{
    print q{Usage: fasta-splitter.pl [options] <file>...
Options:
    --n-parts <N>        - Divide into <N> parts
    --part-size <N>      - Divide into parts of size <N>
    --measure (all|seq|count) - Specify whether all data, sequence length, or
                           number of sequences is used for determining part
                           sizes ('all' by default).
    --line-length        - Set output sequence line length, 0 for single line
                           (default: 60).
    --eol (dos|mac|unix) - Choose end-of-line character ('unix' by default).
    --version            - Show version.
    --help               - Show help.
};
}
