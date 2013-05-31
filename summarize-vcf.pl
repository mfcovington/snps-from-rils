#!/usr/bin/env perl
# summarize-vcf.pl
# Mike Covington
# created: 2013-05-24
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'sum';

my $vcf_file = $ARGV[0];    # "10k.rep_08.E.var.flt.vcf";
open my $vcf_fh, "<", $vcf_file;

my $summary_file = $vcf_file . ".summary";
open my $summary_fh, ">", $summary_file;

while (<$vcf_fh>) {
    next if m|^#|;
    chomp;
    my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, @samples ) =
      split /\t/;
    next unless length ($ref) + length ($alt) == 2;
    my ( $dp4_ref, $dp4_alt ) = $info =~ m/DP4=(\d+,\d+),(\d+,\d+)/;
    dp4_counter($dp4_ref);
    dp4_counter($dp4_alt);
    my $ref_counts = dp4_counter($dp4_ref);
    my $alt_counts = dp4_counter($dp4_alt);
    my $tot_counts = sum $ref_counts, $alt_counts;
    say $summary_fh join "\t", $chr, $pos, $ref, $alt, $ref_counts,
      $alt_counts, $tot_counts;
}

close $vcf_fh;
close $summary_fh;

sub dp4_counter {
    my @counts = split /,/, shift;
    return sum @counts;
}

__END__
A01     84868   .       A       T       96.5    .       DP=13;AF1=0.5368;AC1=90;DP4=3,3,4,3;MQ=30;FQ=99.8;PV4=1,1,0.057,0.00087 GT:PL:GQ        0/1:0,0,0:3     0/1:0,3,37:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:37,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:20,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,3,25:3    0/1:0,0,0:3     0/1:20,3,0:3    0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3     0/1:0,0,0:3
