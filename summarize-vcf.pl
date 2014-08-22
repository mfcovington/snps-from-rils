#!/usr/bin/env perl
# summarize-vcf.pl
# Mike Covington
# created: 2013-05-24
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use List::Util 'sum';
use Getopt::Long;

my $observed_cutoff = 0.1;
my $af1_min         = 0.3;
my $af1_max         = 1 - $af1_min;

my $options = GetOptions (
    "observed_cutoff=f" => \$observed_cutoff,
    "af1_min=f" => \$af1_min,
);

my @vcf_file_list = @ARGV;

for my $vcf_file (@vcf_file_list) {
    open my $vcf_fh, "<", $vcf_file;

    my $summary_file = $vcf_file . ".summary";
    open my $summary_fh, ">", $summary_file;

    my $sample_number;

    while (<$vcf_fh>) {

        # ignore header
        next if m|^##|;

        chomp;
        my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples ) =
          split /\t/;

        # calculate number of samples
        if (m|^#CHROM|) {
            $sample_number = scalar @samples;
            next;
        }

        # ignore INDELs and multiple alternate alleles
        next unless length($ref) + length($alt) == 2;

        my ( $af1, $dp4_ref, $dp4_alt ) =
          $info =~ m/AF1=([^;]+);.+DP4=(\d+,\d+),(\d+,\d+)/;

        # ignore AF1 values too far from 0.5
        next unless $af1 < $af1_max && $af1 > $af1_min;

        my $observed = 0;
        for (@samples) { $observed++ unless m|:0,0,0:|; }
        my $observed_ratio = sprintf "%.2f", $observed / $sample_number;

        # ignore SNPs with coverage in too few samples
        next if $observed_ratio < $observed_cutoff;

        dp4_counter($dp4_ref);
        dp4_counter($dp4_alt);
        my $ref_counts = dp4_counter($dp4_ref);
        my $alt_counts = dp4_counter($dp4_alt);
        my $tot_counts = sum $ref_counts, $alt_counts;
        say $summary_fh join "\t", $chr, $pos, $ref, $alt, $ref_counts,
          $alt_counts, $tot_counts, $af1, $observed_ratio;
    }

    close $vcf_fh;
    close $summary_fh;
}


sub dp4_counter {
    my @counts = split /,/, shift;
    return sum @counts;
}

__END__
A01     84868   .       A       T       96.5    .       DP=13;AF1=0.5368;AC1=90;DP4=3,3,4,3;MQ=30;FQ=99.8;PV4=1,1,0.057,0.00087 GT:PL:GQ        0/1:0,0,0:3     0/1:0,3,37:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:37,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:20,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,3,25:3    0/1:0,0,0:3     0/1:20,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3
