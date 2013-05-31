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

my $count_min;

my %merged;
my %count;
my %conflict;

for my $file ( @vcf_summary_files ) {
    open my $summary_fh, "<", $file;
    # merge_summaries(\$summary_fh, \%merged);
    while ( <$summary_fh> ) {
        my ( $chr, $pos, $ref, $alt ) = split /\t/;
        my $chr_pos = "$chr.$pos";
        my $ref_alt = "$ref.$alt";
        $conflict{$chr_pos} = 1
          if exists $merged{$chr_pos}
          && $merged{$chr_pos} ne $ref_alt;
        $merged{$chr_pos} = $ref_alt;
        $count{$chr_pos}++;
    }
    close $summary_fh;
}

say "merged: ", scalar keys %merged;
say "count: ", scalar keys %count;
say "conflict: ", scalar keys %conflict;

use Data::Printer;

p %conflict;

__END__

$ time ./merge-vcf-summaries.pl A01.rep_01.E.var.flt.vcf.summary A01.rep_02.E.var.flt.vcf.summary A01.rep_03.E.var.flt.vcf.summary
merged: 21174
count: 21174
conflict: 5
{
    A01.11050720   1,
    A01.13174070   1,
    A01.27153976   1,
    A01.27269486   1,
    A01.27741024   1
}

real    0m0.198s
user    0m0.184s
sys 0m0.013s
