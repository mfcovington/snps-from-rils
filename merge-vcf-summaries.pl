#!/usr/bin/env perl
# merge-vcf-summaries.pl
# Mike Covington
# created: 2013-05-30
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'max';

my @vcf_summary_files = @ARGV;
my $replicate_count_min = 2;

my $par1_id = "R500";
my $par2_id = "IMB211";
my $par1_bam = "bwa_tophat_RIL_R500.12-Brapa0830.sorted.dupl_rm.xt_a_u_q20.bam";

#TODO::  incorporate ref-or-alt

my %merged;
my %repeated;
my %conflict;

for my $file (@vcf_summary_files) {
    open my $summary_fh, "<", $file;
    # merge_summaries(\$summary_fh, \%merged);
    while ( <$summary_fh> ) {
        my ( $chr, $pos, $ref, $alt ) = split /\t/;
        my $chr_pos = "$chr.$pos";
        my $ref_alt = "$ref.$alt";
        $conflict{$chr}{$pos} = 1
          if exists $merged{$chr}{$pos}
          && $merged{$chr}{$pos} ne $ref_alt;
        $merged{$chr}{$pos} = $ref_alt;
        $repeated{$chr}{$pos}++;
    }
    close $summary_fh;
}

my %counts;
$counts{merged}   += scalar keys $merged{$_}   for keys %merged;
$counts{repeated} += scalar keys $repeated{$_} for keys %repeated;
$counts{conflict} += scalar keys $conflict{$_} for keys %conflict;
say "merged: ",   $counts{merged};
say "repeated: ", $counts{repeated};
say "conflict: ", $counts{conflict};

for my $chr ( keys %merged ) {
    for my $pos ( keys $merged{$chr}) {
        delete $merged{$chr}{$pos}
          if exists $conflict{$chr}{$pos}
          || $repeated{$chr}{$pos} < $replicate_count_min;
    }
}

$counts{merged_filtered} += scalar keys $merged{$_} for keys %merged;

say "merged (after conflict removal and repeat count filtering): ", $counts{merged_filtered};

use Data::Printer;
p %conflict;

__END__

$ time ./merge-vcf-summaries.pl A01.rep_01.E.var.flt.vcf.summary A01.rep_02.E.var.flt.vcf.summary A01.rep_03.E.var.flt.vcf.summary
merged: 21174
replicate count: 21174
conflict: 5
merged (after conflict removal and replicate count filtering): 12060
{
    A01.11050720   1,
    A01.13174070   1,
    A01.27153976   1,
    A01.27269486   1,
    A01.27741024   1
}

real    0m0.229s
user    0m0.213s
sys 0m0.014s
