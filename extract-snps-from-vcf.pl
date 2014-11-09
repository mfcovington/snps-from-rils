#!/usr/bin/env perl
# summarize-vcf.pl
# Mike Covington
# created: 2013-07-26
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
use Capture::Tiny 'capture_stderr';
use Parallel::ForkManager;


my $threads   = 2;
my $ratio_min = 0.9;    # proportion of reads matching major allele

my $par1_id = "R500";
my $par2_id = "IMB211";

my $sample_dir = "/Users/mfc/git.repos/sample-files/";
my $par1_bam   = "$sample_dir/bam/R500.good.bam";
my $ref_fa     = "$sample_dir/fa/B.rapa_genome_sequence_0830.fa";


my $observed_min  = 0.5;
my $alt_ratio_min = 0.3;
my $alt_ratio_max = 1 - $alt_ratio_min;
my $depth_min     = 4;
my $gq_min        = 13;    # 10 ^ -1.3 == 0.0501
my $het_ratio_max = 0.1;

my $options = GetOptions(
    "observed_min=f" => \$observed_min,
    "alt_ratio_min=f"   => \$alt_ratio_min,
);

my @vcf_file_list = @ARGV;

my $pm_extract_snps = new Parallel::ForkManager($threads);
for my $vcf_file (@vcf_file_list) {
    $pm_extract_snps->start and next;

    open my $vcf_fh, "<", $vcf_file;

    my $summary_file = $vcf_file . ".summary";
    open my $summary_fh, ">", $summary_file;

    my $sample_number;
    my %genotype;

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
        say $summary_fh join "\t", $chr, $pos, $het_count and next if $tot_count == 0;
        # next if $tot_count == 0;
        my $het_ratio = sprintf "%.2f", $het_count / $tot_count;
        say $summary_fh join "\t", $chr, $pos, $het_count, $tot_count, $het_ratio and next unless $het_ratio < $het_ratio_max;
        # next unless $het_ratio < $het_ratio_max;
        my $alt_ratio = sprintf "%.2f", $alt_count / $tot_count;
        next unless $alt_ratio < $alt_ratio_max && $alt_ratio > $alt_ratio_min;
        say $summary_fh join "\t", $chr, $pos, $ref, $alt, $ref_count,
          $alt_count, $tot_count, $af1, $alt_ratio, $observed_ratio, $het_count, $het_ratio;


        $genotype{$chr}{$pos}{ref} = $ref;
        $genotype{$chr}{$pos}{alt} = $alt;

    }

    close $vcf_fh;
    close $summary_fh;

    my $pm_genotype = new Parallel::ForkManager($threads);
    for my $cur_chr ( sort keys %genotype ) {
        $pm_genotype->start and next;

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

            # ignore non-SNP positions
            next unless exists $genotype{$chr}{$pos};

            # ignore insertions
            my @inserts = $bases =~ m|\+|g;
            my $insert_count = scalar @inserts;
            next if $insert_count / $depth > 0.1;

            # get genotype for alternate allele
            my $alt = $genotype{$chr}{$pos}{alt};
            my $alt_genotype =
              get_alt_genotype( $chr, $pos, $ref, $alt, $depth, $bases );

            # ignore ambiguous SNPs & SNPs w/ zero coverage on parent 1
            next if $alt_genotype eq '';

            say $polydb_fh join "\t", $chr, $pos, $ref, $alt, $alt_genotype, 'NA', 'SNP';
        }
        close $mpileup_fh;
        close $polydb_fh;

        $pm_genotype->finish;
    }
    $pm_genotype->wait_all_children;

    $pm_extract_snps->finish;
}
$pm_extract_snps->wait_all_children;


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

