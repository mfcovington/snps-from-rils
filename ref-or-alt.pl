#!/usr/bin/env perl
# ref-or-alt.pl
# Mike Covington
# created: 2013-05-31
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Capture::Tiny 'capture_stderr';

open my $summary_fh, "<", "A01.rep_01.E.var.flt.vcf.summary";
my $count = 0;
while (<$summary_fh>) {
    my ( $chr, $pos, $ref_vcf, $alt ) = split /\t/;
    ref_or_alt( $chr, $pos, $alt );
    exit if $count++ == 10;
}

# my $chr = "A01";
# my $pos = 1537695;
# my $alt = "T";
# my $ref = "C";

sub ref_or_alt {

    my ( $chr, $pos, $alt ) = @_;

    my $par1_id  = "R500";
    my $par1_bam = "R500.good.bam"; # "bwa_tophat_RIL_IMB211.03-Brapa0830.sorted.dupl_rm.xt_a_u_q20.bam";
    my $ref_fa = "B.rapa_genome_sequence_0830.fa";

    my $line;
    capture_stderr {
        $line =
`/Users/mfc/installs/bin/samtools mpileup -r $chr:$pos-$pos -f $ref_fa $par1_bam`;
    };

    return if $line eq "";

    my ( $chr2, $pos2, $ref, $depth, $bases, $quals ) = split /\t/, $line;
    my $skip_count = count_skips($bases);
    my $alt_count  = count_base( $bases, $alt );
    my $ref_count  = count_base($bases);

    print "\n$chr\t$pos\t$ref\t$depth\t$bases\t$quals";
    say "ref count: $ref_count";
    say "alt count: $alt_count";
    say "skip count: $skip_count";
    $depth = $depth - $skip_count;

    if ( $ref_count + $alt_count == 0 ) {
        say "Insufficient coverage at $chr:$pos";
    }
    elsif ( $ref_count == $depth ) {
        say "$par1_id is $ref at $chr:$pos";
    }
    elsif ( $alt_count == $depth ) {
        say "$par1_id is $alt at $chr:$pos";
    }
    else { say "$par1_id is ambiguous at $chr:$pos" }

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

[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
 A01, 8918, C, 15, >>><<>>><<><>aA, IJGGHJJJJIJEH&$

ref count: 13
alt count: 2
R500 is ambiguous at A01:8918
[Finished in 0.4s]
