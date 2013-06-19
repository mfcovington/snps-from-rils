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
use Capture::Tiny 'capture_stderr';
use Parallel::ForkManager;

# TODO:
# - quality score cutoff??

my @vcf_summary_files = @ARGV;

my $verbose             = 0;
my $replicate_count_min = 2;
my $min_ratio           = 0.9;
my $threads             = 3;

my $par1_id  = "R500";
my $par2_id  = "IMB211";
my $par1_bam = "R500.good.bam";
my $ref_fa   = "B.rapa_genome_sequence_0830.fa";

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

for my $chr ( sort keys %merged ) {
    for my $pos ( sort { $a <=> $b } keys $merged{$chr}) {
        delete $merged{$chr}{$pos} and next
          if exists $conflict{$chr}{$pos}
          || $repeated{$chr}{$pos} < $replicate_count_min;
    }
}

$counts{merged_filtered} += scalar keys $merged{$_} for keys %merged;

say "merged (after conflict removal and repeat count filtering): ", $counts{merged_filtered};

use Data::Printer;
p %conflict;



my @chromosomes = qw( A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 );
my $pm = new Parallel::ForkManager($threads);
for my $cur_chr ( @chromosomes ) {
    $pm->start and next;

    my $mpileup_cmd = "/Users/mfc/installs/bin/samtools mpileup -r $cur_chr -f $ref_fa $par1_bam";
    # my $mpileup_cmd = "/Users/mfc/installs/bin/samtools mpileup -r $cur_chr:$pos-$pos -f $ref_fa $par1_bam";
    my $mpileup_fh;
    my $stderr = capture_stderr {    # suppress mpileup output sent to stderr
        open $mpileup_fh,   "-|", $mpileup_cmd;
    };

    # my $mpileup_counter;


    open my $polydb_fh, ">", "polyDB.$cur_chr";
    say $polydb_fh join "\t", 'chr', 'pos', 'ref_base', 'snp_base', 'genotype', 'insert_position', 'SNP_CLASS';

    while ( <$mpileup_fh> ) {
        chomp;
        my ( $chr, $pos, $ref, $depth, $bases, $quals ) = split /\t/;

        next unless exists $merged{$chr}{$pos};
        # $mpileup_counter++;

        my $alt = $merged{$chr}{$pos}{alt};

        my $alt_genotype = get_alt_genotype( $chr, $pos, $ref, $alt, $depth, $bases );
        next if $alt_genotype eq '';
        # say "ALT: $alt_genotype";


        say $polydb_fh join "\t", $chr, $pos, $ref, $alt, $alt_genotype, 'NA', 'SNP';
        # my $nogap_line = <$nogap_fh>;
        # chomp( $gap_line, $nogap_line );
        # my ( $chr, $pos, $gap_cov ) = split /\t/, $gap_line;
        # my $nogap_cov = [ split /\t/, $nogap_line ]->[2];
        # if ( exists $$cov_pos_ref{$chr}{$pos} ) {
        #     $count++;
        #     push @cov_data, [ $sample_id, $chr, $pos, $gap_cov, $nogap_cov ];
        # }

        # populate_and_reset( \$count, \@cov_data, \$schema ) if $count % 100000 == 0;
    }
    close $mpileup_fh;
    close $polydb_fh;

    # say $mpileup_counter;

    $pm->finish;
}
$pm->wait_all_children;


# for my $chr ( sort keys %merged ) {
#     for my $pos ( sort { $a <=> $b } keys $merged{$chr}) {
#     exit;
#     }
# }

        # my $alt = $merged{$chr}{$pos}{alt};
        # my $ref = $merged{$chr}{$pos}{ref};
        # ref_or_alt( $chr, $pos, $alt );
        # say $out_fh join "\t", $chr, $pos, $ref, $alt;

sub get_alt_genotype {

    my ( $chr, $pos, $ref, $alt, $depth, $bases ) = @_;

#     my $par1_id  = "R500";
#     my $par1_bam = "R500.good.bam";
#     my $ref_fa   = "B.rapa_genome_sequence_0830.fa";

#     my $line;
#     capture_stderr {
#         $line =
# `/Users/mfc/installs/bin/samtools mpileup -r $chr:$pos-$pos -f $ref_fa $par1_bam`;
#     };

    # return if $line eq "";

    # my ( $chr2, $pos2, $ref, $depth, $bases, $quals ) = split /\t/, $line;
    my $skip_count  = count_skips($bases);
    my $par1_alt_count   = count_base( $bases, $alt );
    my $par1_ref_count   = count_base($bases);

    if ($verbose) {
        # print "\n$chr\t$pos\t$ref\t$depth\t$bases\t$quals";
        print "\n$chr\t$pos\t$ref\t$depth\t$bases\n";#\t$quals";
        say "ref count:  $par1_ref_count";
        say "alt count:  $par1_alt_count";
        say "skip count: $skip_count";
        $depth = $depth - $skip_count;

        if    ( $depth == 0 )                       { say "Insufficient coverage at $chr:$pos" }
        elsif ( $par1_ref_count >= $min_ratio * $depth ) { say "$par1_id is $ref at $chr:$pos" }
        elsif ( $par1_alt_count >= $min_ratio * $depth ) { say "$par1_id is $alt at $chr:$pos" }
        else                                        { say "$par1_id is ambiguous at $chr:$pos" }
    }

    my $alt_genotype;
    if    ( $depth == 0 )                            { $alt_genotype = '' }
    elsif ( $par1_alt_count >= $min_ratio * $depth ) { $alt_genotype = $par1_id }
    elsif ( $par1_ref_count >= $min_ratio * $depth ) { $alt_genotype = $par2_id }
    else                                             { $alt_genotype = '' }

    return $alt_genotype;
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

