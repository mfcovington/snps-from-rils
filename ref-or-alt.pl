#!/usr/bin/env perl
# ref-or-alt.pl
# Mike Covington
# created: 2013-05-31
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use Capture::Tiny 'capture_stderr';

my $verbose   = 1;
my $min_ratio = 0.9;

open my $summary_fh, "<", "A01.rep_01.E.var.flt.vcf.summary";
my $count = 0;
while (<$summary_fh>) {
    my ( $chr, $pos, $ref_vcf, $alt ) = split /\t/;
    ref_or_alt( $chr, $pos, $alt );
    exit if ++$count == 2;
}

sub ref_or_alt {

    my ( $chr, $pos, $alt ) = @_;

    my $par1_id  = "R500";
    my $par1_bam = "R500.good.bam";
    my $ref_fa   = "B.rapa_genome_sequence_0830.fa";

    my $line;
    capture_stderr {
        $line =
`/Users/mfc/installs/bin/samtools mpileup -r $chr:$pos-$pos -f $ref_fa $par1_bam`;
    };

    return if $line eq "";

    my ( $chr2, $pos2, $ref, $depth, $bases, $quals ) = split /\t/, $line;
    my $skip_count  = count_skips($bases);
    my $alt_count   = count_base( $bases, $alt );
    my $ref_count   = count_base($bases);

    if ($verbose) {
        print "\n$chr\t$pos\t$ref\t$depth\t$bases\t$quals";
        say "ref count:  $ref_count";
        say "alt count:  $alt_count";
        say "skip count: $skip_count";
        $depth = $depth - $skip_count;

        if    ( $depth == 0 )                       { say "Insufficient coverage at $chr:$pos" }
        elsif ( $ref_count >= $min_ratio * $depth ) { say "$par1_id is $ref at $chr:$pos" }
        elsif ( $alt_count >= $min_ratio * $depth ) { say "$par1_id is $alt at $chr:$pos" }
        else                                        { say "$par1_id is ambiguous at $chr:$pos" }
    }
}

sub count_skips {
    my $bases = shift;
    return $bases =~ tr|<>|<>|;
}

sub count_base {
    my ( $bases, $base2count ) = @_;

    my $count;
    if ( defined $base2count ) {
        if    ( $base2count =~ /A/i ) { $count = $bases =~ tr|Aa|Aa| }
        elsif ( $base2count =~ /C/i ) { $count = $bases =~ tr|Cc|Cc| }
        elsif ( $base2count =~ /G/i ) { $count = $bases =~ tr|Gg|Gg| }
        elsif ( $base2count =~ /T/i ) { $count = $bases =~ tr|Tt|Tt| }
    }
    else {
        $count += $bases =~ tr|.,|.,|;
    }

    return $count;
}

__END__

A01 5245    A   4   ,.,.    aaaa
ref count: 4
alt count: 0
skip count: 0
R500 is A at A01:5245

A01 5283    T   5   cCcCC   *.*+!
ref count: 0
alt count: 5
skip count: 0
R500 is C at A01:5283
[Finished in 0.8s]