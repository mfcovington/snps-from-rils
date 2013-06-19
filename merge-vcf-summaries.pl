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
        $conflict{$chr}{$pos} = 1
          if exists $merged{$chr}{$pos}
          && $merged{$chr}{$pos}{alt} ne $alt;
        $merged{$chr}{$pos}{ref} = $ref;
        $merged{$chr}{$pos}{alt} = $alt;
        $repeated{$chr}{$pos}++;
    }
    close $summary_fh;
}

my %counts;
$counts{merged}   += scalar keys $merged{$_}   for keys %merged;
$counts{repeated} += scalar keys $repeated{$_} for keys %repeated;
$counts{conflict} += scalar keys $conflict{$_} for keys %conflict;
say "merged:   ", $counts{merged};
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

$ time ./merge-vcf-summaries.pl A01.rep_*.E.var.flt.vcf.summary
merged:   33149
repeated: 33149
conflict: 13
merged (after conflict removal and repeat count filtering): 23252
{
    A01   {
        5718644    1,
        8030298    1,
        10834898   1,
        11050720   1,
        13174070   1,
        13531866   1,
        21716037   1,
        22878676   1,
        24525675   1,
        27153976   1,
        27154874   1,
        27269486   1,
        27747204   1
    }
}

real    0m0.615s
user    0m0.510s
sys 0m0.026s

