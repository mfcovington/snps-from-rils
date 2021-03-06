#!/usr/bin/env perl
# merge-vcf-summaries.pl
# Mike Covington
# created: 2013-05-30
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use List::Util 'max';
use Capture::Tiny 'capture_stderr';
use Parallel::ForkManager;
use Getopt::Long;

# TODO:
# - quality score cutoff??
# - get chromosome lengths from bam header for use in SNP density stats??
# - usage statement

# RE: Compatibility with VCF files generated
#     by different versions of samtools/bcftools
#
# VCF summaries are DIFFERENT depending on version #!!!
#   - samtools 0.1.18 (r982:295)
#   - samtools 0.1.19-44428cd

my $verbose             = 1;
my $replicate_count_min = 2;      # number of replicates in which SNP is ID'd
my $ratio_min           = 0.9;    # proportion of reads matching major allele
my $threads             = 1;

my ( $par1_id, $par2_id, $par1_bam, $ref_fa, $chr_string, $skip_chr_string );

my $options = GetOptions(
    "par1_id=s"             => \$par1_id,
    "par2_id=s"             => \$par2_id,
    "par1_bam=s"            => \$par1_bam,
    "ref_fa=s"              => \$ref_fa,
    "verbose"               => \$verbose,
    "replicate_count_min=i" => \$replicate_count_min,
    "ratio_min=f"           => \$ratio_min,
    "threads=i"             => \$threads,
    "chr_string=s"          => \$chr_string,
    "skip_chr_string=s"     => \$skip_chr_string,
);

my %chromosomes;
if ( defined $chr_string ) {
    %chromosomes = map { $_ => 1 } split /,/, $chr_string;
}
else {
    %chromosomes = map { $_ => 1 } get_all_chromosomes($par1_bam);
}

if ( defined $skip_chr_string ) {
    delete $chromosomes{$_} for split /,/, $skip_chr_string;
}

my @vcf_summary_files = @ARGV;

my %merged;
my %repeated;
my %conflict;

for my $file (@vcf_summary_files) {
    open my $summary_fh, "<", $file;
    # merge_summaries(\$summary_fh, \%merged);
    while ( <$summary_fh> ) {
        my ( $chr, $pos, $ref, $alt ) = split /\t/;

        next unless exists $chromosomes{$chr};

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
$counts{merged}   += scalar keys %{ $merged{$_} }   for keys %merged;
$counts{repeated} += scalar keys %{ $repeated{$_} } for keys %repeated;

$counts{conflict} = 0;
$counts{conflict} += scalar keys %{ $conflict{$_} } for keys %conflict;

for my $chr ( sort keys %merged ) {
    for my $pos ( sort { $a <=> $b } keys %{ $merged{$chr} } ) {
        delete $merged{$chr}{$pos} and next
          if exists $conflict{$chr}{$pos}
          || $repeated{$chr}{$pos} < $replicate_count_min;
    }
}

$counts{merged_filtered} += scalar keys %{ $merged{$_} } for keys %merged;

if ($verbose) {
    say "merged:   ", $counts{merged};
    say "repeated: ", $counts{repeated};
    say "conflict: ", $counts{conflict};
    say "merged (after conflict removal and repeat count filtering): ",
      $counts{merged_filtered};
}

my $pm = new Parallel::ForkManager($threads);
for my $cur_chr ( sort keys %chromosomes ) {
    $pm->start and next;

    system("samtools index $par1_bam") if ! -e "$par1_bam.bai";
    my $mpileup_cmd = "samtools mpileup -r $cur_chr -f $ref_fa $par1_bam";
    my $mpileup_fh;
    capture_stderr {    # suppress mpileup output sent to stderr
        open $mpileup_fh,   "-|", $mpileup_cmd;
    };

    open my $polydb_fh, ">", "polyDB.$cur_chr";
    say $polydb_fh join "\t", 'chr', 'pos', 'ref_base', 'snp_base', 'genotype', 'insert_position', 'SNP_CLASS';

    while ( <$mpileup_fh> ) {
        chomp;
        my ( $chr, $pos, $ref, $depth, $bases, $quals ) = split /\t/;

        next if $depth == 0;

        # ignore non-SNP positions
        next unless exists $merged{$chr}{$pos};

        # ignore insertions
        my @inserts = $bases =~ m|\+|g;
        my $insert_count = scalar @inserts;
        next if $insert_count / $depth > 0.1;

        # get genotype for alternate allele
        my $alt = $merged{$chr}{$pos}{alt};
        my $alt_genotype =
          get_alt_genotype( $chr, $pos, $ref, $alt, $depth, $bases );

        # ignore ambiguous SNPs & SNPs w/ zero coverage on parent 1
        next if $alt_genotype eq '';

        say $polydb_fh join "\t", $chr, $pos, $ref, $alt, $alt_genotype, 'NA', 'SNP';
    }
    close $mpileup_fh;
    close $polydb_fh;

    $pm->finish;
}
$pm->wait_all_children;

exit;

sub get_alt_genotype {
    my ( $chr, $pos, $ref, $alt, $depth, $bases ) = @_;

    my $skip_count     = count_skips($bases);
    my $par1_alt_count = count_base( $bases, $alt );
    my $par1_ref_count = count_base($bases);

    my $alt_genotype;
    if    ( $depth == 0 )                            { $alt_genotype = '' }
    elsif ( $par1_alt_count >= $ratio_min * $depth ) { $alt_genotype = $par1_id }
    elsif ( $par1_ref_count >= $ratio_min * $depth ) { $alt_genotype = $par2_id }
    else                                             { $alt_genotype = '' }

    return $alt_genotype;
}

sub get_all_chromosomes {
    my $bam = shift;

    my @chromosomes;
    my $samtools_cmd = "samtools view -H $bam";
    open my $bam_header_fh, "-|", $samtools_cmd;
    while (<$bam_header_fh>) {
        next unless /^\@SQ/i;
        /SN:([^\t]+)/;
        push @chromosomes, $1;
    }

    return @chromosomes;
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
