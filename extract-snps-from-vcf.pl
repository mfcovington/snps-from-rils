#!/usr/bin/env perl
# summarize-vcf.pl
# Mike Covington
# created: 2013-07-26
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'sum';
use Getopt::Long;

# my $observed_min  = 0.1;
my $observed_min  = 0.5;
my $alt_ratio_min = 0.3;
my $alt_ratio_max = 1 - $alt_ratio_min;
my $depth_min     = 4;
my $gq_min        = 5;
my $het_ratio_max = 0.1;

my $options = GetOptions(
    "observed_min=f" => \$observed_min,
    "alt_ratio_min=f"   => \$alt_ratio_min,
);

my $vcf_file = $ARGV[0];    # "10k.rep_08.E.var.flt.vcf";
open my $vcf_fh, "<", $vcf_file;

my $summary_file = $vcf_file . ".summary";
open my $summary_fh, ">", $summary_file;

my $sample_number;

while ( my $vcf_line = <$vcf_fh> ) {

    # ignore header
    next if $vcf_line =~ m|^##|;

    # calculate number of samples
    if ($vcf_line =~ m|^#CHROM|) {
        my @col_names = split /\t/, $vcf_line;
        $sample_number =
          scalar @col_names - 9;    # Nine column headers precede sample IDs
        next;
    }

    chomp $vcf_line;
    my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples ) =
      split /\t/, $vcf_line;

    die "Unexpected format column ('$format'). Expected: 'GT:PL:DP:GQ'"
      unless $format eq 'GT:PL:DP:GQ';

    # ignore INDELs and multiple alternate alleles
    next unless length($ref) + length($alt) == 2;

    my ( $af1, $dp4_ref, $dp4_alt ) =
      $info =~ m/AF1=([^;]+);.+DP4=(\d+,\d+),(\d+,\d+)/;

    # ignore AF1 values too far from 0.5
    # next unless $af1 < $alt_ratio_max && $af1 > $alt_ratio_min;

    my $observed_count = 0;
    my $ref_count      = 0;
    my $alt_count      = 0;
    my $het_count      = 0;

    for (@samples) {
        my ( $gt, $pl, $dp, $gq ) = split /:/;
        next if $dp < $depth_min;
        next if $gq < $gq_min;
        # my ( $rr_score, $ra_score, $aa_score ) = split /,/, $pl;
        if    ( $gt eq '0/0' ) { $ref_count++ }
        elsif ( $gt eq '0/1' ) { $het_count++ }    # incorporate???
        elsif ( $gt eq '1/1' ) { $alt_count++ }
        else                   { die "Something wrong with genotype ('$gt')." }
        $observed_count++;
    }




    # for (@samples) { $observed_count++ unless m|:0,0,0:|; }
    my $observed_ratio = sprintf "%.2f", $observed_count / $sample_number;

    # ignore SNPs with coverage in too few samples
    next if $observed_ratio < $observed_min;

    my $tot_count = sum $ref_count, $alt_count;
    # say "oops $pos $het_count" if $tot_count == 0;
    # say $summary_fh join "\t", $chr, $pos, $het_count and next if $tot_count == 0;
    next if $tot_count == 0;
    my $het_ratio = sprintf "%.2f", $het_count / $tot_count;
    # say $summary_fh join "\t", $chr, $pos, $het_count, $tot_count, $het_ratio and next unless $het_ratio < $het_ratio_max;
    next unless $het_ratio < $het_ratio_max;
    my $alt_ratio = sprintf "%.2f", $alt_count / $tot_count;
    next unless $alt_ratio < $alt_ratio_max && $alt_ratio > $alt_ratio_min;
    say $summary_fh join "\t", $chr, $pos, $ref, $alt, $ref_count,
      $alt_count, $tot_count, $af1, $alt_ratio, $observed_ratio, $het_count, $het_ratio;
}

close $vcf_fh;
close $summary_fh;

exit;
